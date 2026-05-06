# Authorship: Gabriel Montenegro de Campos
# Refined by: Gemini
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
import os
from tqdm import tqdm
from sklearn.utils import resample
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import *
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
import xgboost as xgb_lib

# --- Directory Setup ---
outdir = "results_figs"
os.makedirs(outdir, exist_ok=True)

# Load Data
df = pd.read_csv("sheets/Final_df.csv")
colunas_remover = ["Unnamed: 0.1", "Unnamed: 0", "Accession", "Query_ID"] 
df = df.drop(columns=[col for col in colunas_remover if col in df.columns], errors='ignore').fillna(0)

# Data Splitting
X = df.drop('Family', axis=1)
y = df['Family']
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, test_size=0.25, random_state=42)

# --- Bootstrap Analysis Function ---
def run_bootstrap_analysis(model_type="rf", n_bootstrap=1000):
    metrics = {"acc": [], "f1": [], "prec": [], "rec": [], "mcc": []}
    
    label_encoder = LabelEncoder()
    y_train_enc = label_encoder.fit_transform(y_train)
    y_test_enc = label_encoder.transform(y_test)
    
    print(f"\nStarting Bootstrap for {model_type.upper()}...")
    
    for i in tqdm(range(n_bootstrap)):
        X_res, y_res = resample(X_train, y_train_enc, random_state=42 + i)
        
        if model_type == "xgb":
            model = XGBClassifier(eval_metric='logloss', random_state=42+i)
        else:
            model = RandomForestClassifier(class_weight="balanced", random_state=42+i)
            
        model.fit(X_res, y_res)
        y_pred = model.predict(X_test)
        
        metrics["acc"].append(accuracy_score(y_test_enc, y_pred))
        metrics["f1"].append(f1_score(y_test_enc, y_pred, average="weighted"))
        metrics["prec"].append(precision_score(y_test_enc, y_pred, average="weighted", zero_division=0))
        metrics["rec"].append(recall_score(y_test_enc, y_pred, average="weighted"))
        metrics["mcc"].append(matthews_corrcoef(y_test_enc, y_pred))
        
    return metrics

# Run analyses
results_xgb = run_bootstrap_analysis(model_type="xgb")
results_rf = run_bootstrap_analysis(model_type="rf")

# --- Visual Comparison ---
def prepare_plot_data(res_rf, res_xgb):
    data_list = []
    metric_map = {"acc": "Accuracy", "f1": "F1-Score", "prec": "Precision", "rec": "Recall", "mcc": "MCC"}
    
    for m_key, m_name in metric_map.items():
        # Add RF
        temp_rf = pd.DataFrame({'Value': res_rf[m_key]}).assign(Metric=m_name, Classifier='Random Forest')
        # Add XGB
        temp_xgb = pd.DataFrame({'Value': res_xgb[m_key]}).assign(Metric=m_name, Classifier='XGBoost')
        data_list.extend([temp_rf, temp_xgb])
        
    return pd.concat(data_list)

mdf = prepare_plot_data(results_rf, results_xgb)
mdf.to_csv("mdf_metrica_results.csv", index=False)

# Final Boxplot
plt.figure(figsize=(12, 7))
sns.boxplot(x="Metric", y="Value", hue="Classifier", data=mdf, palette="Set2", showmeans=True)
plt.title("Model Performance Comparison (1000 Bootstrap Iterations)")
plt.savefig(f"{outdir}/Model_Comparison_Boxplot.png", dpi=300)
print(f"✅ Comparison plot saved in {outdir}")