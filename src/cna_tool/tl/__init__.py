from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from ..cna_inference import CNAInferer
from ..infer import infer_cnas_from_scrna

def evaluate_cna_call(adata, truth_col='simulated_cnvs', pred_col='cna_profile'):
    #y_true = adata.obs[truth_col].apply(lambda x: 0 if x == '' else 1)
    y_true = adata.obs['simulated_cnvs'].apply(lambda x: 0 if x.strip() == '' else 1)

    y_pred = adata.obs[pred_col].apply(lambda x: 0 if x == '' else 1)

    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    return {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1}


def run_cna_evaluation(adata, params):
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