from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from ..cna_inference import CNAInferer


def evaluate_cna_call(adata, truth_col='simulated_cnvs', pred_col='cna_profile'):
    y_true = adata.obs[truth_col].apply(lambda x: 0 if x=='' else 1)
    y_pred = adata.obs[pred_col].apply(lambda x: 0 if x=='' else 1)

    acc = accuracy_score(y_true, y_pred)
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    return {'accuracy': acc, 'precision': prec, 'recall': rec, 'f1': f1}

def run_cna_evaluation(adata, control, params):
    inferer = CNAInferer(
        adata         = adata.copy(),
        control_adata = control.copy(),
        gtf_df        = None,
        window        = params['window'],
        gain_thr      = params['gain_thr'],
        loss_thr      = params['loss_thr'],
        norm_method   = params['norm_method']
    )
    adata2 = inferer.infer()
    metrics = evaluate_cna_call(adata2)
    return metrics