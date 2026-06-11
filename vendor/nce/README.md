# Introduction 
Elucidation of new compounds based on (mainly) MS2 spectra found by the discovery program in the NewCompoundDiscovery repo, i.e., finding the chemical structure.

# Install and run:
1. Clone the repo.
2. Dependencies are described in the usual files.  The hard one is RDKit.
3. To run the program, you need pickle files of discoveries.  By default it looks for these in
	`c:\cdtemp\<SSID>\<library>\discoveries_<date>.pkl`, and finds the most recent set.  (SSID and library on command line.)
	If the study was an aggregation of sample sets,  e.g.,  {12345, 21345, 11213}, then <SSID> will be formed by sorting the SSIDS and joining with '+',
	e.g., 11213+12345+21345
4. To run programs, start a command prompt, and `cd ./elucidation`.
5. You need ML model files for various parts of the analysis.  Two ways to get these:
	- use `table_builders/synctables.py`. It will get the newest ones from the L: drive.
     It will prompt for every file.  Notice that it will copy IN BOTH DIRECTIONS, if the local files are newer, 
	 so BE CAREFUL.  On the other hand, it never deletes files.
	- rerun the model training programs, which are in `./elucidation/table_builders/\*.py`.
	 Of course, if you do this, you may not have the same models as everyone else.
6. The main program is `elucidator.py`.
	- You can get the command-line options in the usual way:
		`python elucidator.py -h`
	- A typical invocation is:
		`python elucidator.py -s 21676 -l 209`
	- If the study has merged sample sets as above, the invocation is:
		`python elucidator.py -s 11213+12345+21345 -l 209`
	- Individual features:
		`python elucidator.py -s 21676 -l 209 -f F1234 -f F4321`
	- List of individual features in a file, one per line:
		`python elucidator.py -s 21676 -l 209 -F text_file_listing_features`
		
	
7. If you generate some reports for Platform R&D to look at, you need to move them to the `L:` drive.  The `deliver.py` script should no longer be used.
Just copy the directory with the outputs to `L:\New Compound Discovery\<SSID>`.  (The spreadsheets no longer need to be reformatted for the L: drive. Two versions are generated directly.)

8. X-compounds:  The main program is `X_elucidator.py`.
	- You can get the command-line options in the usual way:
		`python X_elucidator.py -h`
	- To elucidate all X compounds:
		`python X_elucidator.py -l 209`
	- To elucidate a single unknown (or a few), the invocation is:
		`python X_elucidator.py -l 209 -c 47942
	- To elucidate a list of comp ids in a file, one per line:
		`python X_elucidator.py -l 209 -C text_file_listing_comp_ids

# Test
Tests are in ./pytest.
