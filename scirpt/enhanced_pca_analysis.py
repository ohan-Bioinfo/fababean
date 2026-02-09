import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
from sklearn.mixture import GaussianMixture
import warnings
warnings.filterwarnings('ignore')
import os

# Set style for professional plots
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 11

def load_passport_data(passport_file):
    """Load and clean passport data"""
    print("Loading passport data...")
    try:
        passport = pd.read_csv(passport_file, sep='\t', encoding='utf-8')
        passport.columns = passport.columns.str.strip()
        print(f"Passport columns: {list(passport.columns)}")
        return passport
    except Exception as e:
        print(f"Error loading passport data: {e}")
        return None

def load_pca_data(eigenvec_file, eigenval_file, passport_data=None):
    """Load PCA results and merge with passport data"""
    eigenvec = pd.read_csv(eigenvec_file, sep='\s+', header=None)
    pc_columns = [f'PC{i+1}' for i in range(len(eigenvec.columns)-2)]
    eigenvec.columns = ['FID', 'IID'] + pc_columns
    
    eigenval = pd.read_csv(eigenval_file, sep='\s+', header=None, names=['Eigenvalue'])
    
    if passport_data is not None:
        print("Merging PCA data with passport information...")
        passport_data['Sent_Code_Clean'] = passport_data['Sent Code'].astype(str).str.strip()
        eigenvec['IID_Clean'] = eigenvec['IID'].astype(str).str.strip()
        
        eigenvec = eigenvec.merge(
            passport_data, 
            left_on='IID_Clean', 
            right_on='Sent_Code_Clean', 
            how='left'
        )
        print(f"Successfully merged {eigenvec['Sent_Code_Clean'].notna().sum()} samples with passport data")
    
    return eigenvec, eigenval

def calculate_variance_explained(eigenval):
    """Calculate variance explained by each principal component"""
    total_variance = eigenval['Eigenvalue'].sum()
    eigenval['Variance_Explained'] = (eigenval['Eigenvalue'] / total_variance) * 100
    eigenval['Cumulative_Variance'] = eigenval['Variance_Explained'].cumsum()
    eigenval['PC'] = [f'PC{i+1}' for i in range(len(eigenval))]
    return eigenval

def find_optimal_clusters(eigenvec, max_k=8):
    """Find optimal number of clusters using multiple metrics"""
    X = eigenvec[['PC1', 'PC2']].values
    n_samples = len(X)
    
    # Adjust max_k based on sample size
    max_k = min(max_k, n_samples - 1)
    
    print(f"Testing cluster numbers from 2 to {max_k}...")
    
    k_range = range(2, max_k + 1)
    results = {
        'k': [],
        'inertia': [],
        'silhouette': [],
        'calinski_harabasz': [],
        'davies_bouldin': [],
        'bic': [],  # For GMM
        'n_clusters': []
    }
    
    # K-means evaluation
    for k in k_range:
        # K-means
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
        kmeans_labels = kmeans.fit_predict(X)
        
        # Calculate metrics
        if k < n_samples:  # Silhouette requires k < n_samples
            silhouette_avg = silhouette_score(X, kmeans_labels)
        else:
            silhouette_avg = np.nan
            
        calinski = calinski_harabasz_score(X, kmeans_labels)
        davies = davies_bouldin_score(X, kmeans_labels)
        
        results['k'].append(k)
        results['inertia'].append(kmeans.inertia_)
        results['silhouette'].append(silhouette_avg)
        results['calinski_harabasz'].append(calinski)
        results['davies_bouldin'].append(davies)
        results['n_clusters'].append(k)
    
    # GMM evaluation
    for k in k_range:
        gmm = GaussianMixture(n_components=k, random_state=42)
        gmm.fit(X)
        results['bic'].append(gmm.bic(X))
    
    results_df = pd.DataFrame(results)
    return results_df, X

def plot_cluster_metrics(results_df, output_dir):
    """Plot multiple cluster evaluation metrics"""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Elbow curve (Inertia)
    ax1.plot(results_df['k'], results_df['inertia'], 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Inertia (Within-cluster sum of squares)', fontsize=12, fontweight='bold')
    ax1.set_title('Elbow Method for Optimal k', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xticks(results_df['k'])
    
    # Annotate the elbow points
    for i, (k, inertia) in enumerate(zip(results_df['k'], results_df['inertia'])):
        ax1.annotate(f'{inertia:.0f}', (k, inertia), 
                    xytext=(0, 10), textcoords='offset points',
                    ha='center', va='bottom', fontsize=9)
    
    # Plot 2: Silhouette Score
    ax2.plot(results_df['k'], results_df['silhouette'], 'ro-', linewidth=2, markersize=8)
    ax2.set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Silhouette Score', fontsize=12, fontweight='bold')
    ax2.set_title('Silhouette Analysis for Optimal k\n(Higher is better)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xticks(results_df['k'])
    
    # Highlight best silhouette score
    best_silhouette_idx = results_df['silhouette'].idxmax()
    best_k_silhouette = results_df.loc[best_silhouette_idx, 'k']
    best_silhouette = results_df.loc[best_silhouette_idx, 'silhouette']
    ax2.plot(best_k_silhouette, best_silhouette, 'g*', markersize=15, 
             label=f'Best: k={best_k_silhouette}, score={best_silhouette:.3f}')
    ax2.legend()
    
    # Plot 3: Calinski-Harabasz Score
    ax3.plot(results_df['k'], results_df['calinski_harabasz'], 'go-', linewidth=2, markersize=8)
    ax3.set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Calinski-Harabasz Score', fontsize=12, fontweight='bold')
    ax3.set_title('Calinski-Harabasz Index\n(Higher is better)', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_xticks(results_df['k'])
    
    # Highlight best Calinski-Harabasz score
    best_calinski_idx = results_df['calinski_harabasz'].idxmax()
    best_k_calinski = results_df.loc[best_calinski_idx, 'k']
    best_calinski = results_df.loc[best_calinski_idx, 'calinski_harabasz']
    ax3.plot(best_k_calinski, best_calinski, 'r*', markersize=15, 
             label=f'Best: k={best_k_calinski}, score={best_calinski:.1f}')
    ax3.legend()
    
    # Plot 4: BIC Score (GMM)
    ax4.plot(results_df['k'], results_df['bic'], 'mo-', linewidth=2, markersize=8)
    ax4.set_xlabel('Number of Clusters (k)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Bayesian Information Criterion (BIC)', fontsize=12, fontweight='bold')
    ax4.set_title('BIC for Gaussian Mixture Models\n(Lower is better)', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xticks(results_df['k'])
    
    # Highlight best BIC score
    best_bic_idx = results_df['bic'].idxmin()
    best_k_bic = results_df.loc[best_bic_idx, 'k']
    best_bic = results_df.loc[best_bic_idx, 'bic']
    ax4.plot(best_k_bic, best_bic, 'b*', markersize=15, 
             label=f'Best: k={best_k_bic}, score={best_bic:.1f}')
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Cluster_Optimization_Metrics.png'), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    return best_k_silhouette, best_k_calinski, best_k_bic

def calculate_elbow_point(inertias, k_values):
    """Calculate the elbow point using the second derivative method"""
    # Calculate first differences
    first_diff = np.diff(inertias)
    # Calculate second differences
    second_diff = np.diff(first_diff)
    # Find the point of maximum curvature (elbow)
    if len(second_diff) > 0:
        elbow_idx = np.argmax(np.abs(second_diff)) + 1  # +1 because we lost one element in diff
        return k_values[elbow_idx]
    return k_values[1]  # Default to 2 if can't determine

def create_detailed_cluster_report(results_df, best_k_silhouette, best_k_calinski, best_k_bic, output_dir):
    """Create a detailed report with cluster recommendations"""
    
    # Calculate elbow point
    elbow_k = calculate_elbow_point(results_df['inertia'].values, results_df['k'].values)
    
    # Create recommendations
    recommendations = {
        'Method': ['Elbow Method', 'Silhouette Score', 'Calinski-Harabasz', 'BIC (GMM)'],
        'Recommended k': [elbow_k, best_k_silhouette, best_k_calinski, best_k_bic],
        'Metric Value': [
            f"{results_df[results_df['k'] == elbow_k]['inertia'].values[0]:.0f}",
            f"{results_df[results_df['k'] == best_k_silhouette]['silhouette'].values[0]:.3f}",
            f"{results_df[results_df['k'] == best_k_calinski]['calinski_harabasz'].values[0]:.1f}",
            f"{results_df[results_df['k'] == best_k_bic]['bic'].values[0]:.1f}"
        ],
        'Interpretation': [
            'Point of maximum curvature in inertia',
            'Highest average silhouette width',
            'Highest variance ratio between clusters',
            'Lowest Bayesian Information Criterion'
        ]
    }
    
    recommendations_df = pd.DataFrame(recommendations)
    
    # Find consensus (most frequent recommended k)
    recommended_ks = [elbow_k, best_k_silhouette, best_k_calinski, best_k_bic]
    consensus_k = max(set(recommended_ks), key=recommended_ks.count)
    
    # Save detailed report
    with open(os.path.join(output_dir, 'Cluster_Analysis_Report.txt'), 'w') as f:
        f.write("FABA BEAN PCA - CLUSTER ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("CLUSTER METRICS SUMMARY:\n")
        f.write("-" * 40 + "\n")
        f.write(results_df.round(4).to_string(index=False))
        f.write("\n\n")
        
        f.write("RECOMMENDED CLUSTER NUMBERS:\n")
        f.write("-" * 40 + "\n")
        f.write(recommendations_df.to_string(index=False))
        f.write("\n\n")
        
        f.write("CONSENSUS RECOMMENDATION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Based on multiple metrics, the recommended number of clusters is: {consensus_k}\n\n")
        
        f.write("INTERPRETATION OF METRICS:\n")
        f.write("-" * 40 + "\n")
        f.write("1. Inertia (Elbow Method): Measures within-cluster sum of squares\n")
        f.write("   - Look for the 'elbow' where inertia stops decreasing rapidly\n")
        f.write("2. Silhouette Score: Measures how similar objects are to their own cluster\n")
        f.write("   - Range: -1 to 1 (higher is better)\n")
        f.write("   > 0.7: Strong structure   0.5-0.7: Reasonable structure\n")
        f.write("   0.25-0.5: Weak structure  < 0.25: No substantial structure\n")
        f.write("3. Calinski-Harabasz: Ratio of between-cluster to within-cluster dispersion\n")
        f.write("   - Higher values indicate better defined clusters\n")
        f.write("4. BIC: Bayesian Information Criterion for Gaussian Mixture Models\n")
        f.write("   - Lower values indicate better model fit\n\n")
        
        f.write("FINAL RECOMMENDATION:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Use k = {consensus_k} for your analysis.\n")
        f.write("This represents the most consistent recommendation across multiple metrics.\n")
    
    # Save metrics to CSV
    results_df.to_csv(os.path.join(output_dir, 'Cluster_Metrics_Detailed.csv'), index=False)
    recommendations_df.to_csv(os.path.join(output_dir, 'Cluster_Recommendations.csv'), index=False)
    
    return consensus_k, recommendations_df

def visualize_all_cluster_solutions(eigenvec, X, max_k, output_dir):
    """Visualize PCA with different cluster numbers for comparison"""
    n_cols = min(3, max_k - 1)
    n_rows = (max_k - 1 + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    for idx, k in enumerate(range(2, max_k + 1)):
        if idx < len(axes):
            ax = axes[idx]
            
            # Perform clustering
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=20)
            labels = kmeans.fit_predict(X)
            
            # Plot
            scatter = ax.scatter(X[:, 0], X[:, 1], c=labels, 
                                cmap='tab10', s=80, alpha=0.8,
                                edgecolor='white', linewidth=1)
            
            # Calculate and display silhouette score
            if k < len(X):
                silhouette_avg = silhouette_score(X, labels)
                ax.set_title(f'k = {k}\nSilhouette: {silhouette_avg:.3f}', 
                           fontsize=12, fontweight='bold')
            else:
                ax.set_title(f'k = {k}', fontsize=12, fontweight='bold')
            
            ax.set_xlabel('PC1')
            ax.set_ylabel('PC2')
            ax.grid(True, alpha=0.3)
    
    # Remove empty subplots
    for idx in range(len(axes)):
        if idx >= (max_k - 1):
            fig.delaxes(axes[idx])
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'All_Cluster_Solutions.png'), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

def apply_final_clustering(eigenvec, optimal_k, output_dir):
    """Apply the final clustering with the optimal k"""
    X = eigenvec[['PC1', 'PC2']].values
    
    # Apply K-means with optimal k
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=20)
    final_labels = kmeans.fit_predict(X)
    eigenvec['Cluster'] = final_labels
    
    # Calculate final metrics
    silhouette_avg = silhouette_score(X, final_labels)
    calinski = calinski_harabasz_score(X, final_labels)
    davies = davies_bouldin_score(X, final_labels)
    
    # Create final cluster assignment file
    cluster_assignments = eigenvec[['IID', 'PC1', 'PC2', 'Cluster']].copy()
    
    # Add passport data if available
    if 'Accession Bank ID' in eigenvec.columns:
        cluster_assignments['Accession_Bank_ID'] = eigenvec['Accession Bank ID']
    if 'Local Names' in eigenvec.columns:
        cluster_assignments['Local_Names'] = eigenvec['Local Names']
    if 'Collection Country' in eigenvec.columns:
        cluster_assignments['Collection_Country'] = eigenvec['Collection Country']
    
    # Save cluster assignments
    cluster_assignments.to_csv(os.path.join(output_dir, f'Final_Cluster_Assignments_k{optimal_k}.csv'), index=False)
    
    # Save cluster statistics
    cluster_stats = cluster_assignments.groupby('Cluster').agg({
        'IID': 'count',
        'PC1': ['mean', 'std'],
        'PC2': ['mean', 'std']
    }).round(4)
    
    cluster_stats.columns = ['Sample_Count', 'PC1_Mean', 'PC1_Std', 'PC2_Mean', 'PC2_Std']
    cluster_stats['Silhouette_Score'] = silhouette_avg
    cluster_stats['Calinski_Harabasz'] = calinski
    cluster_stats['Davies_Bouldin'] = davies
    
    cluster_stats.to_csv(os.path.join(output_dir, f'Final_Cluster_Statistics_k{optimal_k}.csv'))
    
    print(f"\nFinal clustering applied with k={optimal_k}")
    print(f"Silhouette Score: {silhouette_avg:.3f}")
    print(f"Calinski-Harabasz Score: {calinski:.1f}")
    print(f"Davies-Bouldin Score: {davies:.3f}")
    
    return eigenvec

def main():
    """Main function for comprehensive cluster analysis"""
    
    # Set paths
    base_dir = "."
    pca_dir = os.path.join(base_dir, "Diversity/PCA/")
    passport_file = os.path.join(base_dir, "passport.csv")
    
    eigenvec_file = os.path.join(pca_dir, "Faba_PCA.eigenvec")
    eigenval_file = os.path.join(pca_dir, "Faba_PCA.eigenval")
    
    if not os.path.exists(eigenvec_file):
        print(f"Error: {eigenvec_file} not found!")
        return
    
    # Load data
    passport_data = load_passport_data(passport_file) if os.path.exists(passport_file) else None
    eigenvec, eigenval = load_pca_data(eigenvec_file, eigenval_file, passport_data)
    eigenval = calculate_variance_explained(eigenval)
    
    print(f"\nStarting comprehensive cluster analysis...")
    print(f"Total samples: {len(eigenvec)}")
    print(f"PC1 variance: {eigenval.iloc[0]['Variance_Explained']:.2f}%")
    print(f"PC2 variance: {eigenval.iloc[1]['Variance_Explained']:.2f}%")
    
    # Find optimal clusters
    max_k = min(8, len(eigenvec) - 1)  # Maximum k to test
    results_df, X = find_optimal_clusters(eigenvec, max_k)
    
    # Plot metrics
    best_k_silhouette, best_k_calinski, best_k_bic = plot_cluster_metrics(results_df, pca_dir)
    
    # Create detailed report
    consensus_k, recommendations_df = create_detailed_cluster_report(
        results_df, best_k_silhouette, best_k_calinski, best_k_bic, pca_dir
    )
    
    # Visualize all cluster solutions
    visualize_all_cluster_solutions(eigenvec, X, max_k, pca_dir)
    
    # Apply final clustering
    eigenvec = apply_final_clustering(eigenvec, consensus_k, pca_dir)
    
    # Print summary
    print(f"\n{'='*70}")
    print("CLUSTER ANALYSIS COMPLETE - SUMMARY")
    print(f"{'='*70}")
    print(f"Recommended number of clusters: {consensus_k}")
    print(f"\nMetrics for different k values:")
    for k in results_df['k']:
        row = results_df[results_df['k'] == k].iloc[0]
        print(f"k={k}: Inertia={row['inertia']:.0f}, "
              f"Silhouette={row['silhouette']:.3f}, "
              f"Calinski-Harabasz={row['calinski_harabasz']:.1f}")
    
    print(f"\nFiles created in {pca_dir}:")
    print("  - Cluster_Optimization_Metrics.png (visualization of all metrics)")
    print("  - All_Cluster_Solutions.png (PCA plots for all k values)")
    print("  - Cluster_Analysis_Report.txt (detailed analysis and recommendations)")
    print("  - Cluster_Metrics_Detailed.csv (raw metrics for all k)")
    print("  - Cluster_Recommendations.csv (summary of recommendations)")
    print(f"  - Final_Cluster_Assignments_k{consensus_k}.csv (sample cluster assignments)")
    print(f"  - Final_Cluster_Statistics_k{consensus_k}.csv (cluster statistics)")
    
    print(f"\nBased on the analysis, use k={consensus_k} for your clustering.")
    print("You can now provide this cluster number for further visualization.")

if __name__ == "__main__":
    main()