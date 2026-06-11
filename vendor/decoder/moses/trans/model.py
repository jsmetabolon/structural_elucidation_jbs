import torch
import torch.nn as nn
import torch.nn.functional as F


class TranslationModel(nn.Module):
    def __init__(self, vocab, config):
        super().__init__()

        self.vocabulary = vocab
        if hasattr(config, "use_tanh"):
            self.use_tanh = config.use_tanh
        else:
            self.use_tanh = False

        if hasattr(config, "ignore_vae"):
            self.ignore_vae = config.ignore_vae
        else:
            self.ignore_vae = False

        if self.use_tanh:
            assert (
                self.ignore_vae
            ), "--ignore_vae should be set along with --use_tanh"

        # Special symbols
        for ss in ("bos", "eos", "unk", "pad"):
            setattr(self, ss, getattr(vocab, ss))

        # Word embeddings layer
        n_vocab, d_emb = len(vocab), vocab.vectors.size(1)
        self.x_emb = nn.Embedding(n_vocab, d_emb, self.pad)
        self.x_emb.weight.data.copy_(vocab.vectors)
        if config.freeze_embeddings:
            self.x_emb.weight.requires_grad = False

        # Encoder
        if config.q_cell == "gru":
            self.encoder_rnn = nn.GRU(
                d_emb,
                config.q_d_h,
                num_layers=config.q_n_layers,
                batch_first=True,
                dropout=config.q_dropout if config.q_n_layers > 1 else 0,
                bidirectional=config.q_bidir,
            )
        else:
            raise ValueError(
                "Invalid q_cell type, should be one of the ('gru',)"
            )

        q_d_last = config.q_d_h * (2 if config.q_bidir else 1)
        self.q_mu = nn.Linear(q_d_last, config.d_z)
        self.q_logvar = nn.Linear(q_d_last, config.d_z)

        # Decoder
        if config.d_cell == "gru":
            self.decoder_rnn = nn.GRU(
                d_emb + config.d_z,
                config.d_d_h,
                num_layers=config.d_n_layers,
                batch_first=True,
                dropout=config.d_dropout if config.d_n_layers > 1 else 0,
            )
        else:
            raise ValueError(
                "Invalid d_cell type, should be one of the ('gru',)"
            )

        self.decoder_lat = nn.Linear(config.d_z, config.d_d_h)
        self.decoder_fc = nn.Linear(config.d_d_h, n_vocab)
        self.mass_mlp = nn.Sequential(
            nn.Linear(1, config.d_z),
            nn.Tanh(),
            nn.Linear(config.d_z, config.d_z),
        )
        self.ri_mlp = nn.Sequential(
            nn.Linear(1, config.d_z),
            nn.Tanh(),
            nn.Linear(config.d_z, config.d_z),
        )
        self.mass_alpha = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))
        self.ri_alpha = nn.Parameter(torch.tensor(0.0, dtype=torch.float32))

        # Grouping the model's parameters
        if self.use_tanh:
            self.act = nn.Tanh()
            self.encoder = nn.ModuleList(
                [self.encoder_rnn, self.q_mu, self.act, self.q_logvar]
            )
        else:
            self.encoder = nn.ModuleList(
                [self.encoder_rnn, self.q_mu, self.q_logvar]
            )
        self.decoder = nn.ModuleList(
            [self.decoder_rnn, self.decoder_lat, self.decoder_fc]
        )
        self.vae = nn.ModuleList([self.x_emb, self.encoder, self.decoder])

    @property
    def device(self):
        return next(self.parameters()).device

    def string2tensor(self, string, device="model"):
        ids = self.vocabulary.string2ids(string, add_bos=True, add_eos=True)
        tensor = torch.tensor(
            ids,
            dtype=torch.long,
            device=self.device if device == "model" else device,
        )

        return tensor

    def tensor2string(self, tensor):
        ids = tensor.tolist()
        string = self.vocabulary.ids2string(ids, rem_bos=True, rem_eos=True)

        return string

    def _prepare_condition(self, value, batch_size, device):
        if value is None:
            return None
        if not torch.is_tensor(value):
            value = torch.tensor(value, dtype=torch.float32, device=device)
        else:
            value = value.to(device=device, dtype=torch.float32)
        if value.ndim == 0:
            value = value.reshape(1, 1).repeat(batch_size, 1)
        elif value.ndim == 1:
            value = value.unsqueeze(1)
        elif value.ndim != 2:
            raise ValueError(f"Condition tensor must be rank 0, 1, or 2; got shape {tuple(value.shape)}")
        if value.shape[0] != batch_size:
            if value.shape[0] == 1:
                value = value.repeat(batch_size, 1)
            else:
                raise ValueError(
                    f"Condition batch {value.shape[0]} does not match latent batch {batch_size}"
                )
        return value

    def _apply_conditions(self, z, mass=None, ri=None):
        batch_size = z.shape[0]
        mass_cond = self._prepare_condition(mass, batch_size, z.device)
        ri_cond = self._prepare_condition(ri, batch_size, z.device)
        out = z
        if mass_cond is not None:
            out = out + self.mass_alpha * self.mass_mlp(mass_cond)
        if ri_cond is not None:
            out = out + self.ri_alpha * self.ri_mlp(ri_cond)
        return out

    @staticmethod
    def _length_sort_order(lengths):
        return sorted(range(len(lengths)), key=lambda i: lengths[i], reverse=True)

    @staticmethod
    def _is_identity_order(order):
        return all(idx == value for idx, value in enumerate(order))

    @staticmethod
    def _invert_order(order):
        inverse = [0] * len(order)
        for sorted_idx, original_idx in enumerate(order):
            inverse[original_idx] = sorted_idx
        return inverse

    def _restore_original_batch_order(self, tensor, order):
        if self._is_identity_order(order):
            return tensor
        inverse = torch.tensor(
            self._invert_order(order),
            dtype=torch.long,
            device=tensor.device,
        )
        return tensor.index_select(0, inverse)

    def _reconstruction_loss(self, logits, targets):
        if logits.device.type == "mps":
            # `cross_entropy` on MPS is unstable here; move just the loss
            # computation to CPU while keeping gradients connected.
            return F.cross_entropy(
                logits.to("cpu"),
                targets.to("cpu"),
                ignore_index=self.pad,
            )
        return F.cross_entropy(logits, targets, ignore_index=self.pad)

    def forward(self, randomized_smiles, canonical_smiles, mass=None, ri=None):
        """Do the VAE forward step

        :param x: list of tensors of longs, input sentence x
        :return: float, kl term component of loss
        :return: float, recon component of loss
        """
        # Encoder: x -> z, kl_loss
        mu, logvar, z, kl_loss = self.forward_encoder(randomized_smiles)

        # Decoder: x, z -> recon_loss
        if self.ignore_vae:
            # Don't use the sampled latent variable, use the mean
            recon_loss = self.forward_decoder(canonical_smiles, mu, mass=mass, ri=ri)
        else:
            # z is the sampled latent variable from the prior
            recon_loss = self.forward_decoder(canonical_smiles, z, mass=mass, ri=ri)

        # regression_loss = self.regression_loss(mu, logPs)

        return mu, logvar, z, kl_loss, recon_loss

    def forward_encoder(self, x):
        """Encoder step, emulating z ~ E(x) = q_E(z|x)

        :param x: list of tensors of longs, input sentence x
        :return: (n_batch, d_z) of floats, sample of latent vector z
        :return: float, kl term component of loss
        """

        x = [self.x_emb(i_x) for i_x in x]
        lengths = [i_x.size(0) for i_x in x]
        order = self._length_sort_order(lengths)
        if not self._is_identity_order(order):
            x = [x[i] for i in order]
        x = nn.utils.rnn.pack_sequence(x, enforce_sorted=True)

        _, h = self.encoder_rnn(x, None)

        h = h[-(1 + int(self.encoder_rnn.bidirectional)) :]
        h = torch.cat(h.split(1), dim=-1).squeeze(0)
        h = self._restore_original_batch_order(h, order)

        mu, logvar = self.q_mu(h), self.q_logvar(h)

        if self.use_tanh:
            mu = self.act(mu)

        if self.ignore_vae:
            z = mu
            kl_loss = 0
            # logvar = kl_loss = None
        else:
            eps = torch.randn_like(mu)
            z = mu + (logvar / 2).exp() * eps
            kl_loss = 0.5 * (logvar.exp() + mu ** 2 - 1 - logvar).sum(1).mean()
        return mu, logvar, z, kl_loss


    def forward_decoder(self, x, z, mass=None, ri=None):
        """Decoder step, emulating x ~ G(z)

        :param x: list of tensors of longs, input sentence x
        :param z: (n_batch, d_z) of floats, latent vector z
        :return: float, recon component of loss
        """

        z = self._apply_conditions(z, mass=mass, ri=ri)
        lengths = [len(i_x) for i_x in x]
        order = self._length_sort_order(lengths)

        x = nn.utils.rnn.pad_sequence(
            x, batch_first=True, padding_value=self.pad
        )
        x_emb = self.x_emb(x)

        if not self._is_identity_order(order):
            order_index = torch.tensor(order, dtype=torch.long, device=x.device)
            x = x.index_select(0, order_index)
            x_emb = x_emb.index_select(0, order_index)
            z = z.index_select(0, order_index)
            lengths = [lengths[i] for i in order]

        z_0 = z.unsqueeze(1).repeat(1, x_emb.size(1), 1)
        x_input = torch.cat([x_emb, z_0], dim=-1)
        x_input = nn.utils.rnn.pack_padded_sequence(
            x_input, lengths, batch_first=True, enforce_sorted=True
        )

        h_0 = self.decoder_lat(z)
        h_0 = h_0.unsqueeze(0).repeat(self.decoder_rnn.num_layers, 1, 1)

        output, _ = self.decoder_rnn(x_input, h_0)

        output, _ = nn.utils.rnn.pad_packed_sequence(output, batch_first=True)
        y = self.decoder_fc(output)

        recon_loss = self._reconstruction_loss(
            y[:, :-1].contiguous().view(-1, y.size(-1)),
            x[:, 1:].contiguous().view(-1),
        )

        return recon_loss

    def sample_z_prior(self, n_batch):
        """Sampling z ~ p(z) = N(0, I)

        :param n_batch: number of batches
        :return: (n_batch, d_z) of floats, sample of latent z
        """

        return torch.randn(
            n_batch, self.q_mu.out_features, device=self.x_emb.weight.device
        )

    def sample(
        self,
        n_batch,
        max_len=100,
        z=None,
        temp=1.0,
        top_p=0,
        return_log_probs=False,
        mass=None,
        ri=None,
    ):
        """Generating n_batch samples in eval mode (`z` could be
        not on same device)

        :param n_batch: number of sentences to generate
        :param max_len: max len of samples
        :param z: (n_batch, d_z) of floats, latent vector z or None
        :param temp: temperature of softmax
        :return: list of tensors of strings, samples sequence x
        """
        with torch.no_grad():
            if z is None:
                z = self.sample_z_prior(n_batch)
            z = z.to(self.device)
            z = self._apply_conditions(z, mass=mass, ri=ri)
            z_0 = z.unsqueeze(1)

            # Initial values
            h = self.decoder_lat(z)
            h = h.unsqueeze(0).repeat(self.decoder_rnn.num_layers, 1, 1)
            w = torch.tensor(self.bos, device=self.device).repeat(n_batch)
            x = torch.tensor([self.pad], device=self.device).repeat(
                n_batch, max_len
            )
            x[:, 0] = self.bos
            end_pads = torch.tensor([max_len], device=self.device).repeat(
                n_batch
            )
            eos_mask = torch.zeros(n_batch, dtype=torch.bool, device=self.device)
            log_prob_sums = torch.zeros(
                n_batch, dtype=torch.float32, device=self.device
            )
            token_counts = torch.zeros(
                n_batch, dtype=torch.float32, device=self.device
            )

            # Generating cycle
            for i in range(1, max_len):
                x_emb = self.x_emb(w).unsqueeze(1)
                x_input = torch.cat([x_emb, z_0], dim=-1)

                o, h = self.decoder_rnn(x_input, h)
                y = self.decoder_fc(o.squeeze(1))
                # TODO  y = top_k_top_p_filtering(y, top_p=top_p)
                y = F.softmax(y / temp, dim=-1)

                w = torch.multinomial(y, 1)[:, 0]
                chosen = torch.gather(y, 1, w.unsqueeze(1)).squeeze(1).clamp_min(1e-12)
                active_mask = (~eos_mask).to(dtype=torch.float32)
                log_prob_sums = log_prob_sums + torch.log(chosen) * active_mask
                token_counts = token_counts + active_mask
                x[~eos_mask, i] = w[~eos_mask]
                # i_eos_mask = ~eos_mask & (w == self.eos)
                i_eos_mask = (~eos_mask) & (w == self.eos)
                end_pads[i_eos_mask] = i + 1
                eos_mask = eos_mask | i_eos_mask

            # Converting `x` to list of tensors
            new_x = []
            for i in range(x.size(0)):
                new_x.append(x[i, : end_pads[i]])
            decoded = [self.tensor2string(i_x) for i_x in new_x]
            if return_log_probs:
                mean_log_probs = (
                    log_prob_sums / token_counts.clamp_min(1.0)
                ).detach().cpu().tolist()
                return decoded, mean_log_probs
            return decoded
