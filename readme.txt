
The script for running predictions using eSIG-Net is as follows:
bash test.sh


test.sh is an example prediction script, h5file.py is the feature generation script and model_pred_opt_fold_0.pth is the pre-trained model checkpoint file.

Briefly, h5file.py converts custom FASTA sequences into the required SDNN .h5 feature file. Each FASTA record ID will be used as the key in the .h5 file, so the IDs in the input prediction CSV should match the FASTA IDs exactly.

The prediction CSV should follow the same six-column format without a header:

mutant_protein_id, partner_protein_id, label_after, wildtype_protein_id, label_before, mutation

For example:

RAD51D.68.GtoK,XRCC2,0,RAD51D,1,68G>K
RAD51D.239.RtoP,XRCC2,0,RAD51D,1,239R>P
RAD51D.68.GtoW,XRCC2,0,RAD51D,1,68G>W

If users want to perform prediction but do not have wildtype reference PPI labels, label_before (wildtype PPI) can be derived based on experimental evidence (e.g., human interactome databases) or computational prediction (e.g., AlphaFold-Multimer).
