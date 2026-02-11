import pandas as pd
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, accuracy_score

# --- CONFIGURATIONS ---
DATA_PATH = "data/Final_df.csv"
MODEL_OUTPUT = "rf_model.pkl"
RANDOM_STATE = 42

def load_and_preprocess_data(file_path):
    """
    Loads the dataset and performs initial cleaning.
    """
    print(f"Loading data from: {file_path}")
    df = pd.read_csv(file_path)
    
    # Columns to drop (metadata or identifiers)
    columns_to_remove = [
        "Unnamed: 0.1", "Unnamed: 0", "Accession", "Organism_Name", 
        "Species", "Genus", "Molecule_type", "Country", "Host", "Query_ID"
    ]
    
    # Drop only if columns exist in the dataframe
    df = df.drop(columns=[col for col in columns_to_remove if col in df.columns])
    
    # Handle missing values
    df = df.fillna(0)
    
    return df

def train_and_evaluate(df, target='Family'):
    """
    Splits the data, trains the Random Forest model, and prints metrics.
    """
    X = df.drop(target, axis=1)
    y = df[target]

    # Stratified split to maintain class proportions
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, stratify=y, test_size=0.25, random_state=RANDOM_STATE
    )

    print(f"Training set size: {X_train.shape[0]} samples")
    
    # Initialize and fit the model
    rf_clf = RandomForestClassifier(
        class_weight="balanced", 
        random_state=RANDOM_STATE, 
        n_jobs=-1
    )
    rf_clf.fit(X_train, y_train)

    # Predictions and Metrics
    y_pred = rf_clf.predict(X_test)
    
    print("\n--- Model Evaluation ---")
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))

    return rf_clf

if __name__ == "__main__":
    try:
        # Step 1: Data Preparation
        clean_df = load_and_preprocess_data(DATA_PATH)
        
        # Step 2: Training
        model = train_and_evaluate(clean_df)
        
        # Step 3: Saving the model
        joblib.dump(model, MODEL_OUTPUT)
        print(f"\n? Model successfully saved to {MODEL_OUTPUT}")
        
    except FileNotFoundError:
        print(f"? Error: The file '{DATA_PATH}' was not found. Please check the path.")
    except Exception as e:
        print(f"? An unexpected error occurred: {e}")