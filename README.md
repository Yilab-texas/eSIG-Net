# eSIG-Net

Accurate prediction of single-mutation induced perturbations on protein interactions using a language model.

## Overview

eSIG-Net (edgetic mutation Sequence-based Interaction Grammar Network) is a novel sequence-based "Interaction Language Model" for predicting protein interaction alterations caused by single mutations. eSIG-Net combines various protein sequence embeddings, introduces a mutation-encoding module with syntax and evolutionary insights, and employs contrastive learning to evaluate mutation-induced interaction changes. eSIG-Net significantly outperforms current state-of-the-art sequence-based and structure-based prediction methods at predicting mutational impact on protein interactions.

## Installation

### Requirements

- Python 3.7+
- PyTorch
- CUDA (for GPU acceleration, recommended)

### Dependencies

Install required packages:

```bash
pip install torch torchvision torchaudio
pip install numpy pandas scikit-learn
pip install h5py pyyaml optuna
```

Or create a conda environment:

```bash
conda create -n esignet python=3.8
conda activate esignet
conda install pytorch torchvision torchaudio -c pytorch
pip install numpy pandas scikit-learn h5py pyyaml optuna
```

## Data Setup

### Download Datasets

The datasets for eSIG-Net model training and evaluation can be downloaded from:
**[Zenodo Repository](https://zenodo.org/records/15486737)**

### Data Structure

After downloading and extracting the datasets, organize your data as follows:

```
datasets/
├── csvs/
│   ├── germline.csv
│   ├── random_fold_0.csv
│   ├── random_fold_1.csv
│   ├── random_fold_2.csv
│   ├── random_fold_3.csv
│   └── random_fold_4.csv
└── embeddings/
    ├── sdnn_corrected_ppi.h5
    └── 650M_1/
        ├── smile.10015.pt
        ├── smile.10015_A309T.pt
        ├── smile.10015_G429S.pt
        └── ...
```

### Data Format

- **CSV files**: Each CSV file should contain columns (no header) for:
  - `id_A`: Original protein A identifier
  - `id_B`: Protein B identifier
  - `label`: Interaction label for original pair
  - `id_A_2`: Mutated protein A identifier
  - `label_2`: Interaction label for mutated pair
  - `mute_idx`: Mutation index/position

- **H5 file**: SDNN embeddings in HDF5 format with protein IDs as keys

- **PT files**: ESM embeddings stored as PyTorch tensors with structure:
  - `['representations'][33]`: ESM layer 33 embeddings

## Configuration

The training configuration is specified in `config_train.yaml`. Key parameters include:

- `train_task`: Model task type (default: `sdnn`)
- `exec_mode`: Execution mode (`train` or `eval`)
- `k_fold`: Number of folds for cross-validation (default: 5)
- `batch_size_trn`: Training batch size (default: 32)
- `n_epochs`: Number of training epochs (default: 8)
- `lr_init`: Initial learning rate (default: 0.001)
- `dropout`: Dropout probability (default: 0.5)
- Model architecture parameters (`fc_hidden`, `pred_layer_num`, etc.)

## Usage

### Training

#### Basic Training (Single Fold)

Train a single fold:

```bash
python main.py --config_fname config_train.yaml --fold 0 --save_model --save_results --mdl_dpath tmp_results/exp01
```

#### K-Fold Cross-Validation

Train all folds (specified in config):

```bash
python main.py --config_fname config_train.yaml --save_model --save_results --mdl_dpath tmp_results/exp01
```

The model will automatically train on all folds specified by `k_fold` in the config file.

### Evaluation

Evaluate a trained model:

```bash
python main.py --config_fname config_train.yaml \
  --eval \
  --eval_csv datasets/ppi_case/case_reports_mutation.csv \
  --eval_model_fpath tmp_results/ODrd01_5fold/exp0/model_pred_opt_fold_0.pth \
  --feat_h5_fpath datasets/embeddings/sdnn/sdnn_corrected_ppi.h5 \
  --save_result_fpath eval_results/ODrd01_5fold_exp0_fold0.csv \
  --save_feature_fpath eval_results/feat_ODrd01_5fold_exp0_fold0.pt \
  --fold 0
```

### Using the Training Script

Alternatively, use the provided shell script:

```bash
bash train.sh
```

Note: The script includes GLIBC path configurations that may need adjustment for your system.

## Output

### Training Outputs

- **Model checkpoints**: Saved in `mdl_dpath` directory as `model_pred_opt_fold_{fold}.pth`
- **Results CSV**: Contains training metrics and predictions (if `save_results=True`)

### Evaluation Outputs

- **Results CSV**: Contains predictions and metrics for each sample
- **Feature files**: Extracted features saved as PyTorch tensors (.pt format)

## Project Structure

```
eSIG-Net/
├── main.py                 # Main entry point
├── config_train.yaml       # Training configuration
├── ppi_learner.py          # Model training and evaluation logic
├── ppi_dataset.py          # Dataset loading and preprocessing
├── losses.py               # Loss functions
├── metric_recorder.py      # Metrics tracking
├── backbones/
│   └── sdnn/
│       └── sdnn_model.py   # SDNN model architecture
├── utils/
│   └── sdnn_utils.py       # Utility functions
└── datasets/               # Data directory (to be created)
    ├── csvs/
    └── embeddings/
```

## Notes

- Ensure CUDA is available if using GPU acceleration
- The model requires both SDNN embeddings (H5 format) and ESM embeddings (PT format)
- For k-fold cross-validation, ensure all fold CSV files (`random_fold_0.csv` through `random_fold_{k-1}.csv`) are present
- ESM embeddings should be stored in `datasets/embeddings/650M_1/` directory

## Citation

If you use eSIG-Net in your research, please cite:

```
eSIG-Net: Accurate prediction of single-mutation induced perturbations on
protein interactions using a language model
Xingxin Pan, Aditya Shrawat, Sidharth Raghavan, Yuntao Yang, Zhao Li,
W. Jim Zheng, S. Gail Eckhardt, Erxi Wu, Juan I. Fuxman Bass, Daniel F. Jarosz,
Sidi Chen, Jason H. Huang, Daniel J. McGrail, Nidhi Sahni, and S. Stephen Yi
```

**Correspondence:** Xingxin Pan (xingxin.pan@bswhealth.org), Nidhi Sahni (nidhi.sahni.2025@gmail.com), Daniel J. McGrail (mcgraid@ccf.org), or S. Stephen Yi (song.yi@bcm.edu)
