#some other possible ligands not pursued for now.
obabel -:"CC(C(=O)O)SSC(C(=O)O)C" -opdb --gen3D > 2-MPA-SS.pdb
obabel -:"CC(C)([C@H](C(=O)O)N)SSC(C)([C@H](C(=O)O)N)C" -opdb --gen3D > Pen-SS.pdb
obabel -:"C(C(=O)O)SSC(C(=O)O)" -opdb --gen3D > MCA-ss.pdb
obabel -:"C(C(=O)O)CCCCSSCCCCC(C(=O)O)" -opdb --gen3D > 6-MHA-ss.pdb
obabel -:"C1=CC=C(C(=C1)C(=O)O)SSC(C(=C2)C(=O)O)=CC=C2" -opdb --gen3D > TSA-ss.pdb
obabel -:"C1=CC(=CC=C1C(=O)O)SSC(=CC=C2C(=O)O)C=C2" -opdb --gen3D > 4-MBA-ss.pdb
obabel -:"[H][C@](N)(CSSC[C@]([H])(N)C(O)=O)C(O)=O" -opdb --gen3D > cystine.pdb
