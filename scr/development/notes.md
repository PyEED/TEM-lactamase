

### AlphaFoold SetUp Test

`docker run --rm --runtime=nvidia --gpus all nvidia/cuda:11.8.0-base-ubuntu20.04 nvidia-smi`
(alphafold_env) nab@harry:~/Niklas$ python alphafold/docker/run_docker.py     --fasta_paths=test.fasta     --max_template_date=2022-01-01     --data_dir=/media/database/alphafold     --output_dir=/home/nab/Niklas
I0208 13:55:14.536303 124881755502400 run_docker.py:116] Mounting /home/nab/Niklas -> /mnt/fasta_path_0
I0208 13:55:14.536455 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/uniref90 -> /mnt/uniref90_database_path
I0208 13:55:14.536557 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/mgnify -> /mnt/mgnify_database_path
I0208 13:55:14.536630 124881755502400 run_docker.py:116] Mounting /media/database/alphafold -> /mnt/data_dir
I0208 13:55:14.536696 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/pdb_mmcif/mmcif_files -> /mnt/template_mmcif_dir
I0208 13:55:14.536764 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/pdb_mmcif -> /mnt/obsolete_pdbs_path
I0208 13:55:14.536912 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/pdb70 -> /mnt/pdb70_database_path
I0208 13:55:14.536986 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/uniref30 -> /mnt/uniref30_database_path
I0208 13:55:14.537058 124881755502400 run_docker.py:116] Mounting /media/database/alphafold/bfd -> /mnt/bfd_database_path
I0208 13:55:15.161474 124881755502400 run_docker.py:259] /bin/bash: /opt/conda/lib/libtinfo.so.6: no version information available (required by /bin/bash)
I0208 13:55:18.921923 124881755502400 run_docker.py:259] I0208 13:55:18.921152 138467799745152 templates.py:858] Using precomputed obsolete pdbs /mnt/obsolete_pdbs_path/obsolete.dat.
I0208 13:55:19.765541 124881755502400 run_docker.py:259] I0208 13:55:19.764944 138467799745152 xla_bridge.py:863] Unable to initialize backend 'rocm': NOT_FOUND: Could not find registered platform with name: "rocm". Available platform names are: CUDA
I0208 13:55:19.766190 124881755502400 run_docker.py:259] I0208 13:55:19.765954 138467799745152 xla_bridge.py:863] Unable to initialize backend 'tpu': INTERNAL: Failed to open libtpu.so: libtpu.so: cannot open shared object file: No such file or directory
I0208 13:55:25.150957 124881755502400 run_docker.py:259] I0208 13:55:25.150285 138467799745152 run_alphafold.py:524] Have 5 models: ['model_1_pred_0', 'model_2_pred_0', 'model_3_pred_0', 'model_4_pred_0', 'model_5_pred_0']
I0208 13:55:25.151198 124881755502400 run_docker.py:259] I0208 13:55:25.150457 138467799745152 run_alphafold.py:538] Using random seed 106187891190118892 for the data pipeline
I0208 13:55:25.151337 124881755502400 run_docker.py:259] I0208 13:55:25.150712 138467799745152 run_alphafold.py:245] Predicting test
I0208 13:55:25.151792 124881755502400 run_docker.py:259] I0208 13:55:25.151488 138467799745152 jackhmmer.py:133] Launching subprocess "/usr/bin/jackhmmer -o /dev/null -A /tmp/tmpku59z6ig/output.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 /mnt/fasta_path_0/test.fasta /mnt/uniref90_database_path/uniref90.fasta"
I0208 13:55:25.152986 124881755502400 run_docker.py:259] I0208 13:55:25.152644 138467799745152 utils.py:36] Started Jackhmmer (uniref90.fasta) query
I0208 14:03:03.086571 124881755502400 run_docker.py:259] I0208 14:03:03.085747 138467799745152 utils.py:40] Finished Jackhmmer (uniref90.fasta) query in 457.933 seconds
I0208 14:03:03.309976 124881755502400 run_docker.py:259] I0208 14:03:03.309392 138467799745152 jackhmmer.py:133] Launching subprocess "/usr/bin/jackhmmer -o /dev/null -A /tmp/tmpvp5vm1qt/output.sto --noali --F1 0.0005 --F2 5e-05 --F3 5e-07 --incE 0.0001 -E 0.0001 --cpu 8 -N 1 /mnt/fasta_path_0/test.fasta /mnt/mgnify_database_path/mgy_clusters_2022_05.fa"
I0208 14:03:03.310821 124881755502400 run_docker.py:259] I0208 14:03:03.310449 138467799745152 utils.py:36] Started Jackhmmer (mgy_clusters_2022_05.fa) query
I0208 14:15:25.082044 124881755502400 run_docker.py:259] I0208 14:15:25.081380 138467799745152 utils.py:40] Finished Jackhmmer (mgy_clusters_2022_05.fa) query in 741.771 seconds
I0208 14:15:27.065240 124881755502400 run_docker.py:259] I0208 14:15:27.064711 138467799745152 hhsearch.py:85] Launching subprocess "/usr/bin/hhsearch -i /tmp/tmp0dot88ld/query.a3m -o /tmp/tmp0dot88ld/output.hhr -maxseq 1000000 -d /mnt/pdb70_database_path/pdb70"
I0208 14:15:27.066130 124881755502400 run_docker.py:259] I0208 14:15:27.065907 138467799745152 utils.py:36] Started HHsearch query
I0208 14:15:56.405305 124881755502400 run_docker.py:259] I0208 14:15:56.404733 138467799745152 utils.py:40] Finished HHsearch query in 29.339 seconds
I0208 14:15:57.167103 124881755502400 run_docker.py:259] I0208 14:15:57.166532 138467799745152 hhblits.py:128] Launching subprocess "/usr/bin/hhblits -i /mnt/fasta_path_0/test.fasta -cpu 4 -oa3m /tmp/tmp1hl3ktjw/output.a3m -o /dev/null -n 3 -e 0.001 -maxseq 1000000 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -d /mnt/bfd_database_path/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt -d /mnt/uniref30_database_path/UniRef30_2021_03"
I0208 14:15:57.167811 124881755502400 run_docker.py:259] I0208 14:15:57.167536 138467799745152 utils.py:36] Started HHblits query
I0208 14:24:38.790810 124881755502400 run_docker.py:259] I0208 14:24:38.747937 138467799745152 utils.py:40] Finished HHblits query in 521.575 seconds
I0208 14:24:38.935492 124881755502400 run_docker.py:259] I0208 14:24:38.933544 138467799745152 templates.py:879] Searching for template for: MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW
I0208 14:24:39.144154 124881755502400 run_docker.py:259] I0208 14:24:39.142736 138467799745152 templates.py:267] Found an exact template match 4zj1_A.
I0208 14:24:39.604608 124881755502400 run_docker.py:259] I0208 14:24:39.603835 138467799745152 templates.py:267] Found an exact template match 2p74_B.
I0208 14:24:41.200925 124881755502400 run_docker.py:259] I0208 14:24:41.200066 138467799745152 templates.py:267] Found an exact template match 4bd0_A.
I0208 14:24:41.357366 124881755502400 run_docker.py:259] I0208 14:24:41.356810 138467799745152 templates.py:267] Found an exact template match 6nfd_A.
I0208 14:24:41.659289 124881755502400 run_docker.py:259] I0208 14:24:41.658745 138467799745152 templates.py:267] Found an exact template match 1m40_A.
I0208 14:24:41.817106 124881755502400 run_docker.py:259] I0208 14:24:41.816545 138467799745152 templates.py:267] Found an exact template match 1n9b_A.
I0208 14:24:42.048118 124881755502400 run_docker.py:259] I0208 14:24:42.047436 138467799745152 templates.py:267] Found an exact template match 2b5r_B.
I0208 14:24:42.466316 124881755502400 run_docker.py:259] I0208 14:24:42.465743 138467799745152 templates.py:267] Found an exact template match 4ua6_A.
I0208 14:24:42.555942 124881755502400 run_docker.py:259] I0208 14:24:42.555326 138467799745152 templates.py:267] Found an exact template match 1g6a_A.
I0208 14:24:42.818588 124881755502400 run_docker.py:259] I0208 14:24:42.817541 138467799745152 templates.py:267] Found an exact template match 6afo_B.
I0208 14:24:43.260471 124881755502400 run_docker.py:259] I0208 14:24:43.259819 138467799745152 templates.py:267] Found an exact template match 6td0_A.
I0208 14:24:43.431005 124881755502400 run_docker.py:259] I0208 14:24:43.430433 138467799745152 templates.py:267] Found an exact template match 1o7e_B.
I0208 14:24:43.668414 124881755502400 run_docker.py:259] I0208 14:24:43.667691 138467799745152 templates.py:267] Found an exact template match 6niq_A.
I0208 14:24:43.825849 124881755502400 run_docker.py:259] I0208 14:24:43.824931 138467799745152 templates.py:267] Found an exact template match 4mbh_A.
I0208 14:24:44.194959 124881755502400 run_docker.py:259] I0208 14:24:44.194534 138467799745152 templates.py:267] Found an exact template match 5ne2_B.
I0208 14:24:44.419879 124881755502400 run_docker.py:259] I0208 14:24:44.419458 138467799745152 templates.py:267] Found an exact template match 6qwb_A.
I0208 14:24:44.579257 124881755502400 run_docker.py:259] I0208 14:24:44.578722 138467799745152 templates.py:267] Found an exact template match 6c7a_A.
I0208 14:24:44.813825 124881755502400 run_docker.py:259] I0208 14:24:44.813270 138467799745152 templates.py:286] Found a fuzzy sequence-only match 6dmh_A.
I0208 14:24:45.286695 124881755502400 run_docker.py:259] I0208 14:24:45.286113 138467799745152 templates.py:267] Found an exact template match 4c75_D.
I0208 14:24:45.415336 124881755502400 run_docker.py:259] I0208 14:24:45.414770 138467799745152 templates.py:267] Found an exact template match 6bn3_A.
I0208 14:24:46.167034 124881755502400 run_docker.py:259] I0208 14:24:46.166486 138467799745152 pipeline.py:234] Uniref90 MSA size: 10000 sequences.
I0208 14:24:46.167183 124881755502400 run_docker.py:259] I0208 14:24:46.166624 138467799745152 pipeline.py:235] BFD MSA size: 2460 sequences.
I0208 14:24:46.167282 124881755502400 run_docker.py:259] I0208 14:24:46.166654 138467799745152 pipeline.py:236] MGnify MSA size: 501 sequences.
I0208 14:24:46.167401 124881755502400 run_docker.py:259] I0208 14:24:46.166685 138467799745152 pipeline.py:237] Final (deduplicated) MSA size: 12900 sequences.
I0208 14:24:46.167842 124881755502400 run_docker.py:259] I0208 14:24:46.166904 138467799745152 pipeline.py:239] Total number of templates (NB: this can include bad templates and is later filtered to top 4): 20.
I0208 14:24:46.231397 124881755502400 run_docker.py:259] I0208 14:24:46.230882 138467799745152 run_alphafold.py:276] Running model model_1_pred_0 on test
I0208 14:24:50.838901 124881755502400 run_docker.py:259] I0208 14:24:50.837518 138467799745152 model.py:165] Running predict with shape(feat) = {'aatype': (4, 286), 'residue_index': (4, 286), 'seq_length': (4,), 'template_aatype': (4, 4, 286), 'template_all_atom_masks': (4, 4, 286, 37), 'template_all_atom_positions': (4, 4, 286, 37, 3), 'template_sum_probs': (4, 4, 1), 'is_distillation': (4,), 'seq_mask': (4, 286), 'msa_mask': (4, 508, 286), 'msa_row_mask': (4, 508), 'random_crop_to_size_seed': (4, 2), 'template_mask': (4, 4), 'template_pseudo_beta': (4, 4, 286, 3), 'template_pseudo_beta_mask': (4, 4, 286), 'atom14_atom_exists': (4, 286, 14), 'residx_atom14_to_atom37': (4, 286, 14), 'residx_atom37_to_atom14': (4, 286, 37), 'atom37_atom_exists': (4, 286, 37), 'extra_msa': (4, 5120, 286), 'extra_msa_mask': (4, 5120, 286), 'extra_msa_row_mask': (4, 5120), 'bert_mask': (4, 508, 286), 'true_msa': (4, 508, 286), 'extra_has_deletion': (4, 5120, 286), 'extra_deletion_value': (4, 5120, 286), 'msa_feat': (4, 508, 286, 49), 'target_feat': (4, 286, 22)}
I0208 14:26:57.844202 124881755502400 run_docker.py:259] I0208 14:26:57.843049 138467799745152 model.py:175] Output shape was {'distogram': {'bin_edges': (63,), 'logits': (286, 286, 64)}, 'experimentally_resolved': {'logits': (286, 37)}, 'masked_msa': {'logits': (508, 286, 23)}, 'predicted_lddt': {'logits': (286, 50)}, 'structure_module': {'final_atom_mask': (286, 37), 'final_atom_positions': (286, 37, 3)}, 'plddt': (286,), 'ranking_confidence': ()}
I0208 14:26:57.844417 124881755502400 run_docker.py:259] I0208 14:26:57.843188 138467799745152 run_alphafold.py:288] Total JAX model model_1_pred_0 on test predict time (includes compilation time, see --benchmark): 127.0s
I0208 14:26:58.044636 124881755502400 run_docker.py:259] I0208 14:26:58.044024 138467799745152 run_alphafold.py:276] Running model model_2_pred_0 on test
I0208 14:27:01.635061 124881755502400 run_docker.py:259] I0208 14:27:01.632167 138467799745152 model.py:165] Running predict with shape(feat) = {'aatype': (4, 286), 'residue_index': (4, 286), 'seq_length': (4,), 'template_aatype': (4, 4, 286), 'template_all_atom_masks': (4, 4, 286, 37), 'template_all_atom_positions': (4, 4, 286, 37, 3), 'template_sum_probs': (4, 4, 1), 'is_distillation': (4,), 'seq_mask': (4, 286), 'msa_mask': (4, 508, 286), 'msa_row_mask': (4, 508), 'random_crop_to_size_seed': (4, 2), 'template_mask': (4, 4), 'template_pseudo_beta': (4, 4, 286, 3), 'template_pseudo_beta_mask': (4, 4, 286), 'atom14_atom_exists': (4, 286, 14), 'residx_atom14_to_atom37': (4, 286, 14), 'residx_atom37_to_atom14': (4, 286, 37), 'atom37_atom_exists': (4, 286, 37), 'extra_msa': (4, 1024, 286), 'extra_msa_mask': (4, 1024, 286), 'extra_msa_row_mask': (4, 1024), 'bert_mask': (4, 508, 286), 'true_msa': (4, 508, 286), 'extra_has_deletion': (4, 1024, 286), 'extra_deletion_value': (4, 1024, 286), 'msa_feat': (4, 508, 286, 49), 'target_feat': (4, 286, 22)}
I0208 14:28:34.324337 124881755502400 run_docker.py:259] I0208 14:28:34.323230 138467799745152 model.py:175] Output shape was {'distogram': {'bin_edges': (63,), 'logits': (286, 286, 64)}, 'experimentally_resolved': {'logits': (286, 37)}, 'masked_msa': {'logits': (508, 286, 23)}, 'predicted_lddt': {'logits': (286, 50)}, 'structure_module': {'final_atom_mask': (286, 37), 'final_atom_positions': (286, 37, 3)}, 'plddt': (286,), 'ranking_confidence': ()}
I0208 14:28:34.324619 124881755502400 run_docker.py:259] I0208 14:28:34.323359 138467799745152 run_alphafold.py:288] Total JAX model model_2_pred_0 on test predict time (includes compilation time, see --benchmark): 92.7s
I0208 14:28:34.503692 124881755502400 run_docker.py:259] I0208 14:28:34.503018 138467799745152 run_alphafold.py:276] Running model model_3_pred_0 on test
I0208 14:28:37.690137 124881755502400 run_docker.py:259] I0208 14:28:37.688968 138467799745152 model.py:165] Running predict with shape(feat) = {'aatype': (4, 286), 'residue_index': (4, 286), 'seq_length': (4,), 'is_distillation': (4,), 'seq_mask': (4, 286), 'msa_mask': (4, 512, 286), 'msa_row_mask': (4, 512), 'random_crop_to_size_seed': (4, 2), 'atom14_atom_exists': (4, 286, 14), 'residx_atom14_to_atom37': (4, 286, 14), 'residx_atom37_to_atom14': (4, 286, 37), 'atom37_atom_exists': (4, 286, 37), 'extra_msa': (4, 5120, 286), 'extra_msa_mask': (4, 5120, 286), 'extra_msa_row_mask': (4, 5120), 'bert_mask': (4, 512, 286), 'true_msa': (4, 512, 286), 'extra_has_deletion': (4, 5120, 286), 'extra_deletion_value': (4, 5120, 286), 'msa_feat': (4, 512, 286, 49), 'target_feat': (4, 286, 22)}
I0208 14:29:55.803868 124881755502400 run_docker.py:259] I0208 14:29:55.803133 138467799745152 model.py:175] Output shape was {'distogram': {'bin_edges': (63,), 'logits': (286, 286, 64)}, 'experimentally_resolved': {'logits': (286, 37)}, 'masked_msa': {'logits': (512, 286, 23)}, 'predicted_lddt': {'logits': (286, 50)}, 'structure_module': {'final_atom_mask': (286, 37), 'final_atom_positions': (286, 37, 3)}, 'plddt': (286,), 'ranking_confidence': ()}
I0208 14:29:55.804049 124881755502400 run_docker.py:259] I0208 14:29:55.803262 138467799745152 run_alphafold.py:288] Total JAX model model_3_pred_0 on test predict time (includes compilation time, see --benchmark): 78.1s
I0208 14:29:55.984359 124881755502400 run_docker.py:259] I0208 14:29:55.983707 138467799745152 run_alphafold.py:276] Running model model_4_pred_0 on test
I0208 14:29:59.030320 124881755502400 run_docker.py:259] I0208 14:29:59.029333 138467799745152 model.py:165] Running predict with shape(feat) = {'aatype': (4, 286), 'residue_index': (4, 286), 'seq_length': (4,), 'is_distillation': (4,), 'seq_mask': (4, 286), 'msa_mask': (4, 512, 286), 'msa_row_mask': (4, 512), 'random_crop_to_size_seed': (4, 2), 'atom14_atom_exists': (4, 286, 14), 'residx_atom14_to_atom37': (4, 286, 14), 'residx_atom37_to_atom14': (4, 286, 37), 'atom37_atom_exists': (4, 286, 37), 'extra_msa': (4, 5120, 286), 'extra_msa_mask': (4, 5120, 286), 'extra_msa_row_mask': (4, 5120), 'bert_mask': (4, 512, 286), 'true_msa': (4, 512, 286), 'extra_has_deletion': (4, 5120, 286), 'extra_deletion_value': (4, 5120, 286), 'msa_feat': (4, 512, 286, 49), 'target_feat': (4, 286, 22)}
I0208 14:31:13.643562 124881755502400 run_docker.py:259] I0208 14:31:13.642972 138467799745152 model.py:175] Output shape was {'distogram': {'bin_edges': (63,), 'logits': (286, 286, 64)}, 'experimentally_resolved': {'logits': (286, 37)}, 'masked_msa': {'logits': (512, 286, 23)}, 'predicted_lddt': {'logits': (286, 50)}, 'structure_module': {'final_atom_mask': (286, 37), 'final_atom_positions': (286, 37, 3)}, 'plddt': (286,), 'ranking_confidence': ()}
I0208 14:31:13.643769 124881755502400 run_docker.py:259] I0208 14:31:13.643102 138467799745152 run_alphafold.py:288] Total JAX model model_4_pred_0 on test predict time (includes compilation time, see --benchmark): 74.6s
I0208 14:31:13.822230 124881755502400 run_docker.py:259] I0208 14:31:13.821843 138467799745152 run_alphafold.py:276] Running model model_5_pred_0 on test
I0208 14:31:16.839077 124881755502400 run_docker.py:259] I0208 14:31:16.838348 138467799745152 model.py:165] Running predict with shape(feat) = {'aatype': (4, 286), 'residue_index': (4, 286), 'seq_length': (4,), 'is_distillation': (4,), 'seq_mask': (4, 286), 'msa_mask': (4, 512, 286), 'msa_row_mask': (4, 512), 'random_crop_to_size_seed': (4, 2), 'atom14_atom_exists': (4, 286, 14), 'residx_atom14_to_atom37': (4, 286, 14), 'residx_atom37_to_atom14': (4, 286, 37), 'atom37_atom_exists': (4, 286, 37), 'extra_msa': (4, 1024, 286), 'extra_msa_mask': (4, 1024, 286), 'extra_msa_row_mask': (4, 1024), 'bert_mask': (4, 512, 286), 'true_msa': (4, 512, 286), 'extra_has_deletion': (4, 1024, 286), 'extra_deletion_value': (4, 1024, 286), 'msa_feat': (4, 512, 286, 49), 'target_feat': (4, 286, 22)}
I0208 14:32:30.091397 124881755502400 run_docker.py:259] I0208 14:32:30.090774 138467799745152 model.py:175] Output shape was {'distogram': {'bin_edges': (63,), 'logits': (286, 286, 64)}, 'experimentally_resolved': {'logits': (286, 37)}, 'masked_msa': {'logits': (512, 286, 23)}, 'predicted_lddt': {'logits': (286, 50)}, 'structure_module': {'final_atom_mask': (286, 37), 'final_atom_positions': (286, 37, 3)}, 'plddt': (286,), 'ranking_confidence': ()}
I0208 14:32:30.091645 124881755502400 run_docker.py:259] I0208 14:32:30.090906 138467799745152 run_alphafold.py:288] Total JAX model model_5_pred_0 on test predict time (includes compilation time, see --benchmark): 73.3s
I0208 14:32:35.518680 124881755502400 run_docker.py:259] I0208 14:32:35.517703 138467799745152 amber_minimize.py:178] alterations info: {'nonstandard_residues': [], 'removed_heterogens': set(), 'missing_residues': {}, 'missing_heavy_atoms': {}, 'missing_terminals': {<Residue 285 (TRP) of chain 0>: ['OXT']}, 'Se_in_MET': [], 'removed_chains': {0: []}}
I0208 14:32:35.715488 124881755502400 run_docker.py:259] I0208 14:32:35.714844 138467799745152 amber_minimize.py:408] Minimizing protein, attempt 1 of 100.
I0208 14:32:36.026923 124881755502400 run_docker.py:259] I0208 14:32:36.026144 138467799745152 amber_minimize.py:69] Restraining 2212 / 4439 particles.
I0208 14:32:38.472040 124881755502400 run_docker.py:259] I0208 14:32:38.471008 138467799745152 amber_minimize.py:178] alterations info: {'nonstandard_residues': [], 'removed_heterogens': set(), 'missing_residues': {}, 'missing_heavy_atoms': {}, 'missing_terminals': {}, 'Se_in_MET': [], 'removed_chains': {0: []}}
I0208 14:32:41.428249 124881755502400 run_docker.py:259] I0208 14:32:41.427688 138467799745152 amber_minimize.py:500] Iteration completed: Einit 371101.85 Efinal -7155.45 Time 1.28 s num residue violations 0 num residue exclusions 0
I0208 14:32:42.129069 124881755502400 run_docker.py:259] I0208 14:32:42.128388 138467799745152 run_alphafold.py:414] Final timings for test: {'features': 1761.0301377773285, 'process_features_model_1_pred_0': 4.606117248535156, 'predict_and_compile_model_1_pred_0': 127.00606513023376, 'process_features_model_2_pred_0': 3.5878937244415283, 'predict_and_compile_model_2_pred_0': 92.691326379776, 'process_features_model_3_pred_0': 3.1857049465179443, 'predict_and_compile_model_3_pred_0': 78.11442136764526, 'process_features_model_4_pred_0': 3.04539155960083, 'predict_and_compile_model_4_pred_0': 74.61389017105103, 'process_features_model_5_pred_0': 3.0162715911865234, 'predict_and_compile_model_5_pred_0': 73.25267243385315, 'relax_model_1_pred_0': 11.269187450408936}
(alphafold_env) nab@harry:~/Niklas$ ls
DockerShare  TEM-lactamase  alphafold  non-work  pyeed  share  test  test.fasta
(alphafold_env) nab@harry:~/Niklas$ cd test
(alphafold_env) nab@harry:~/Niklas/test$ ls
confidence_model_1_pred_0.json  features.pkl  ranked_1.pdb  ranked_4.cif                relaxed_model_1_pred_0.pdb  result_model_5_pred_0.pkl     unrelaxed_model_2_pred_0.pdb  unrelaxed_model_5_pred_0.cif
confidence_model_2_pred_0.json  msas          ranked_2.cif  ranked_4.pdb                result_model_1_pred_0.pkl   timings.json                  unrelaxed_model_3_pred_0.cif  unrelaxed_model_5_pred_0.pdb
confidence_model_3_pred_0.json  ranked_0.cif  ranked_2.pdb  ranking_debug.json          result_model_2_pred_0.pkl   unrelaxed_model_1_pred_0.cif  unrelaxed_model_3_pred_0.pdb
confidence_model_4_pred_0.json  ranked_0.pdb  ranked_3.cif  relax_metrics.json          result_model_3_pred_0.pkl   unrelaxed_model_1_pred_0.pdb  unrelaxed_model_4_pred_0.cif
confidence_model_5_pred_0.json  ranked_1.cif  ranked_3.pdb  relaxed_model_1_pred_0.cif  result_model_4_pred_0.pkl   unrelaxed_model_2_pred_0.cif  unrelaxed_model_4_pred_0.pdb
(alphafold_env) nab@harry:~/Niklas/test$ cat confidence_model_1_pred_0.json 
{"residueNumber":[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286],"confidenceScore":[44.09,47.95,47.59,51.35,57.15,56.18,52.49,55.03,57.11,52.9,63.34,65.07,59.09,61.02,58.4,57.84,56.33,54.48,55.41,46.35,51.36,53.13,55.76,56.42,62.32,90.76,95.94,96.53,97.46,98.11,98.48,98.25,98.14,98.52,98.33,98.3,98.36,98.54,98.42,98.46,98.52,98.82,98.86,98.87,98.79,98.77,98.42,98.43,98.38,97.54,97.56,97.46,98.0,98.29,97.98,97.68,98.38,98.69,98.54,98.54,98.65,98.7,98.84,98.91,98.87,98.92,98.86,98.85,98.93,98.92,98.91,98.94,98.93,98.94,98.9,98.87,98.86,98.89,98.9,98.76,98.79,98.82,98.71,98.52,98.49,98.58,98.69,98.77,98.75,98.67,98.79,98.77,98.84,98.72,98.73,98.43,98.24,98.39,98.63,98.74,98.66,97.94,97.81,98.63,98.66,98.86,98.85,98.62,98.73,98.76,98.45,98.24,98.57,98.74,98.89,98.89,98.9,98.9,98.91,98.95,98.92,98.92,98.87,98.92,98.89,98.81,98.81,98.57,98.86,98.86,98.91,98.94,98.93,98.92,98.9,98.94,98.91,98.75,98.26,98.53,98.56,98.74,98.84,98.69,98.81,98.9,98.86,98.77,98.76,98.85,98.62,98.28,97.75,98.31,98.76,98.66,98.84,98.91,98.9,98.92,98.85,98.86,98.84,98.86,98.67,98.69,98.67,98.73,98.53,98.42,98.53,97.98,94.93,98.48,98.63,98.83,98.89,98.95,98.94,98.89,98.88,98.76,98.88,98.94,98.87,98.76,98.88,98.86,98.69,98.53,98.72,98.48,97.63,95.19,96.63,97.02,98.65,98.7,98.51,98.58,98.8,98.68,98.64,98.81,98.87,98.61,98.78,98.86,98.88,98.76,98.7,98.56,98.32,97.64,98.51,98.41,98.56,98.75,98.64,98.58,98.34,97.76,97.91,97.61,96.85,97.35,98.51,98.71,98.8,98.87,98.91,98.91,98.86,98.79,98.75,98.45,98.53,98.52,98.34,98.87,98.87,98.85,98.92,98.94,98.93,98.89,98.82,98.58,98.13,97.18,96.98,97.78,98.25,98.42,98.76,98.86,98.9,98.94,98.93,98.93,98.88,98.76,98.27,98.18,98.1,98.19,98.52,98.42,98.49,98.6,98.76,98.78,98.66,98.66,98.78,98.68,98.44,98.53,98.32,98.06,98.0,97.67,97.08,95.72,94.2,93.69],"confidenceCategory":["D","D","D","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","L","D","L","L","L","L","L","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H"]}(alphafold_env) nab@harry:~/Niklas/test$ 

in der run.py     [docker.types.DeviceRequest(driver="nvidia", capabilities=[["gpu"]], count=-1)]
docker file in zwei zeilen bei conda ainstalltion


# **NFS Setup Guide**

## **Overview**
This document provides step-by-step instructions for setting up an **NFS server** on `129.69.129.131` and mounting it on an **NFS client** at `129.69.129.130`.

### **NFS Server**: `129.69.129.131`
- Shared Directory: `/srv/nfs/shared_folder`

### **NFS Client**: `129.69.129.130`
- Mount Point: `/mnt/nfs_shared`

---

## **Step 1: Configure the NFS Server (`129.69.129.131`)**

### **1. Install the NFS Server Package**
```bash
sudo apt update && sudo apt install nfs-kernel-server -y  # Debian/Ubuntu
sudo yum install nfs-utils -y  # RHEL/CentOS
```

### **2. Create and Configure the Shared Directory**
```bash
sudo mkdir -p /srv/nfs/shared_folder
sudo chown -R nobody:nogroup /srv/nfs/shared_folder
sudo chmod -R 777 /srv/nfs/shared_folder
```

### **3. Configure NFS Exports**
Edit the exports file:
```bash
sudo nano /etc/exports
```
Add the following line:
```
/srv/nfs/shared_folder 129.69.129.130(rw,sync,no_subtree_check)
```

### **4. Apply Changes and Restart the NFS Server**
```bash
sudo exportfs -ra
sudo systemctl restart nfs-server
```

### **5. Allow NFS in the Firewall (If Applicable)**
```bash
sudo ufw allow from 129.69.129.130 to any port nfs  # Ubuntu/Debian
sudo firewall-cmd --permanent --add-service=nfs  # RHEL/CentOS
sudo firewall-cmd --reload
```

---

## **Step 2: Configure the NFS Client (`129.69.129.130`)**

### **1. Install the NFS Client Package**
```bash
sudo apt update && sudo apt install nfs-common -y  # Debian/Ubuntu
sudo yum install nfs-utils -y  # RHEL/CentOS
```

### **2. Create a Mount Point**
```bash
sudo mkdir -p /mnt/nfs_shared
```

### **3. Mount the NFS Share**
```bash
sudo mount -t nfs 129.69.129.131:/srv/nfs/shared_folder /mnt/nfs_shared
```

### **4. Verify the Mount**
```bash
df -h | grep nfs
ls -l /mnt/nfs_shared
```

### **5. Make the Mount Persistent Across Reboots**
Edit the `/etc/fstab` file:
```bash
sudo nano /etc/fstab
```
Add this line:
```
129.69.129.131:/srv/nfs/shared_folder  /mnt/nfs_shared  nfs  rw,sync  0  0
```

Apply the changes:
```bash
sudo mount -a
```

### **6. Test Writing a File**
```bash
touch /mnt/nfs_shared/testfile.txt
echo "NFS works!" > /mnt/nfs_shared/testfile.txt
ls -l /mnt/nfs_shared
```

---

## **Unmounting the NFS Share (If Needed)**
If you need to unmount the NFS share manually:
```bash
sudo umount /mnt/nfs_shared
```
If the mount is busy:
```bash
sudo umount -l /mnt/nfs_shared
```

---

## **Summary**
This guide sets up an NFS share on `129.69.129.131` and allows `129.69.129.130` to mount and access it. The steps ensure proper installation, configuration, and automatic mounting for persistent access.

🚀 **Now, your NFS setup is ready!**

