


python main.py \
  --config_fname config_eval5.yaml \
  --eval \
  --eval_csv datasets/csvs/topredcition.csv \
  --eval_model_fpath ./model_pred_opt_fold_0.pth \
  --feat_h5_fpath datasets/embeddings/combin3_dedup.h5 \
  --save_result_fpath eval_results/topredcition.csv \
  --save_feature_fpath eval_results/topredcition.pt


