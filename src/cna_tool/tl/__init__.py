from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from ..cna_inference import CNAInferer
from ..infer import infer_cnas_from_scrna

def evaluate_cna_call(adata, truth_col='simulated_cnvs', pred_col='cna_profile'):
    """
    Evaluates the performance of CNA prediction by comparing predicted CNA profiles
    with simulated ground truth labels.

    Parameters:
    - adata: AnnData object containing simulated ground truth and predicted CNA profiles.
    - truth_col (str): Column name in adata.obs with the simulated ground truth CNA labels.
                       Empty string means no CNA; non-empty means presence of CNA.
    - pred_col (str): Column name in adata.obs with predicted CNA profile strings.

    Returns:
    - dict: Dictionary containing accuracy, precision, recall, and F1 score.
    """
    #y_true = adata.obs[truth_col].apply(lambda x: 0 if x == '' else 1)
    y_true = adata.obs['simulated_cnvs'].apply(lambda x: 0 if x.strip() == '' else 1)

    y_pred = adata.obs[pred_col].apply(lambda x: 0 if x == '' else 1)

    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    return {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1}


def run_cna_evaluation(adata, params):
    """
    Runs the full CNA inference and evaluation pipeline using the provided parameters.

    Parameters:
    - adata: AnnData object with simulated CNAs (used as ground truth).
    - params (dict): Dictionary of CNA inference parameters including:
        - 'window': Sliding window size
        - 'gain_thr': Gain threshold
        - 'loss_thr': Loss threshold
        - 'norm_method': Normalization method ('log2_ratio' or 'zscore')

    Returns:
    - dict: Evaluation metrics (accuracy, precision, recall, F1) after running CNA inference.
    """
    adata_eval = infer_cnas_from_scrna(
        adata.copy(),
        cna_label_col='simulated_cnvs',
        diploid_labels=[''],
        gtf_df=None,
        window=params['window'],
        gain_thr=params['gain_thr'],
        loss_thr=params['loss_thr'],
        norm_method=params['norm_method']
    )
    return evaluate_cna_call(adata_eval)