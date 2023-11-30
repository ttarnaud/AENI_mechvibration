### Threshold 05 : RMS threshold used in strength curve calculations is 5%
v2: Amp_fkApdpset_NCfixed4_nos4l
version 2 => kAus is fixed to 2 times value of kAus when nr of dipoles highest 10^5

Difference Input_cEG_param_kAus_dps16_nos4l_v2 vs Input_cEG_param_nos4l_v1:
new sim method: calculate field for single vibration amplitude then based on this calculate for others (using information gained on linearity of DataBase study. That result linear is in vibration amplitude for fixed directions)
less dipoles per run, maxtime decreased

sOSC: test if results with single oscilater (only tagged dipole) differ highly from results where the 2 times kAus artificial field was applied. Huge difference in computation time

### Threshold 10: RMS threshold used in strength curve calculations is 10%
From this point on (03/03/2022) threshold is included in Input file

### aIfun: alternative input function. Now train of alpha's except single Alpha,
Load input function containing Taus, Alphadelays, and Amps see aTInput_220308.mat (tis method works best to ensure same current is used when compared to results in calculating RMS)