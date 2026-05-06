import pandas as pd
import joblib
import os

# Path
base_path = "." 
model_path = "models"
data_path = "spades/DarkMatterFiles"
output_path = "results_predictions"

os.makedirs(output_path, exist_ok=True)

#Load models
xgb = joblib.load(os.path.join(model_path, "modelo_xgb.pkl"))
rf = joblib.load(os.path.join(model_path, "modelo_rf.pkl"))
le = joblib.load(os.path.join(model_path, "label_encoder.pkl"))

# Which Sample would you like to explore (example: FN13)
sample_id = "FN13"
df_sample = pd.read_csv(os.path.join(data_path, f"{sample_id}_DarkMatter.csv"))

model_features = xgb.get_booster().feature_names
missing_cols = [col for col in model_features if col not in df_sample.columns]

for col in missing_cols:
    df_sample[col] = 0

X = df_sample[model_features]

pfam_cols = [c for c in X.columns if c.startswith('PF')] 
mask_has_annotation = (X[pfam_cols] != 0).any(axis=1)

X_filtered = X.loc[mask_has_annotation].copy()
nodes_filtered = df_sample.loc[mask_has_annotation, 'Node']

print(f"{sample_id}: {len(X)} contigs -> {len(X_filtered)} anottaded.")

# Predctions
print("Predicting...")
preds_xgb = xgb.predict(X_filtered)
preds_rf = rf.predict(X_filtered)

probs_xgb = xgb.predict_proba(X_filtered)
probs_rf = rf.predict_proba(X_filtered)

resultados = pd.DataFrame({
    'Sample': sample_id,
    'Node': nodes_filtered,
    'Prediction_XGB': le.inverse_transform(preds_xgb),
    'Prediction_RF': le.inverse_transform(preds_rf)
})

resultados['XGB'] = probs_xgb.max(axis=1)
resultados['RF'] = probs_rf.max(axis=1)

resultados.to_csv(os.path.join(output_path, f"{sample_id}_final_predictions.csv"), index=False)
print(f"✅ Resultados salvos em {output_path}")