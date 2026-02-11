import pandas as pd
import joblib
import xgboost as xgb
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, accuracy_score

# --- CONFIGURATIONS ---
DATA_PATH = "data/Final_df.csv"
MODEL_OUTPUT = "xgb_model.pkl"
ENCODER_OUTPUT = "label_encoder.pkl"
RANDOM_STATE = 42

def load_and_preprocess_data(file_path):
    """
    Loads the dataset and removes unnecessary metadata columns.
    """
    print(f"Reading data from: {file_path}")
    df = pd.read_csv(file_path)
    
    # Metadata columns to be removed
    drop_cols = [
        "Unnamed: 0.1", "Unnamed: 0", "Accession", "Organism_Name", 
        "Species", "Genus", "Molecule_type", "Country", "Host", "Query_ID"
    ]
    
    df = df.drop(columns=[c for c in drop_cols if c in df.columns])
    return df.fillna(0)

def train_xgb_model(df, target_col='Family'):
    """
    Encodes labels, splits data, and trains the XGBoost Classifier.
    """
    X = df.drop(target_col, axis=1)
    y = df[target_col]

    # XGBoost requires numerical labels for multiclass classification
    encoder = LabelEncoder()
    y_encoded = encoder.fit_transform(y)
    
    # Save the encoder to decode predictions later
    joblib.dump(encoder, ENCODER_OUTPUT)
    print(f"? Label Encoder saved to {ENCODER_OUTPUT}")

    X_train, X_test, y_train, y_test = train_test_split(
        X, y_encoded, stratify=y_encoded, test_size=0.25, random_state=RANDOM_STATE
    )

    print(f"Training XGBoost version: {xgb.__version__}")
    
    # Initialize and train
    xgb_clf = XGBClassifier(
        eval_metric='logloss', 
        random_state=RANDOM_STATE,
        use_label_encoder=False
    )
    
    xgb_clf.fit(X_train, y_train)

    # Evaluation
    y_pred = xgb_clf.predict(X_test)
    print("\n--- Model Evaluation ---")
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    return xgb_clf

if __name__ == "__main__":
    try:
        # Flow execution
        data = load_and_preprocess_data(DATA_PATH)
        model = train_xgb_model(data)
        
        # Save model
        joblib.dump(model, MODEL_OUTPUT)
        print(f"? XGBoost model saved successfully: {MODEL_OUTPUT}")
        
    except FileNotFoundError:
        print(f"? Error: File not found at {DATA_PATH}")
    except Exception as e:
        print(f"? An error occurred: {e}")