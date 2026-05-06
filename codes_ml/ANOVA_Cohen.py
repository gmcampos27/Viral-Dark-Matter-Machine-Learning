import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.formula.api import ols
import os

# --- Funções Estatísticas ---

def calculate_cohen_d(g1, g2):
    """(Cohen's d)."""
    n1, n2 = len(g1), len(g2)
    var1, var2 = np.var(g1, ddof=1), np.var(g2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    return (np.mean(g1) - np.mean(g2)) / pooled_std if pooled_std > 0 else 0

def interpret_cohen(d):
    d_abs = abs(d)
    if d_abs < 0.2: return "Negligible"
    if d_abs < 0.5: return "Small"
    if d_abs < 0.8: return "Medium"
    return "Large"


input_file = "mdf_metrica_results.csv"
if not os.path.exists(input_file):
    print(f"❌ Error {input_file} not found!")
else:
    mdf = pd.read_csv(input_file)
    
    tabela_resumo = []
    tabela_suplementar = []


    for metrica in mdf["metrica"].unique():
        df_m = mdf[mdf["metrica"] == metrica]
        
        modelo = ols("value ~ C(Classifier)", data=df_m).fit()
        anova_results = sm.stats.anova_lm(modelo, typ=2)
        p_anova = anova_results["PR(>F)"][0]
        
        classifiers = df_m["Classifier"].unique()
        if len(classifiers) == 2:
            g1_name, g2_name = classifiers[0], classifiers[1]
            g1 = df_m[df_m["Classifier"] == g1_name]["value"]
            g2 = df_m[df_m["Classifier"] == g2_name]["value"]
            
            d_val = calculate_cohen_d(g1, g2)
            eff_size = interpret_cohen(d_val)
        else:
            g1_name, g2_name, d_val, eff_size = classifiers[0], "N/A", 0, "N/A"

        tabela_resumo.append({
            "Metric": metrica,
            "Comparison": f"{g1_name} vs {g2_name}",
            "p-value": f"{p_anova:.2e}",
            "Cohen's d": round(d_val, 4),
            "Effect Size": eff_size
        })

        for source in ["C(Classifier)", "Residual"]:
            tabela_suplementar.append({
                "Metric": metrica,
                "Source": source.replace("C(Classifier)", "Model"),
                "SS": round(anova_results.loc[source, "sum_sq"], 6),
                "df": int(anova_results.loc[source, "df"]),
                "F-value": round(anova_results.loc[source, "F"], 4) if source != "Residual" else np.nan,
                "p-value": anova_results.loc[source, "PR(>F)"] if source != "Residual" else np.nan
            })


    df_resumo = pd.DataFrame(tabela_resumo)
    df_resumo.to_csv("Model_Comparison_Statistical_Summary.csv", index=False)

    df_detalhes = pd.DataFrame(tabela_suplementar)
    df_detalhes.to_csv("Model_Comparison_ANOVA_Details.csv", index=False)

    print("\nTables:")
    print("1. Model_Comparison_Statistical_Summary.csv (Resumo p-value e Cohen)")
    print("2. Model_Comparison_ANOVA_Details.csv (Detalhes de soma de quadrados)\n")
    print(df_resumo)