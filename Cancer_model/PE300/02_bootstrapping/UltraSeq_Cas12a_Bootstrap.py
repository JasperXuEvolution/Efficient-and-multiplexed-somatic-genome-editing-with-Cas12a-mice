#!/usr/bin/env python
# coding: utf-8

# # Deep sequencing analysis
# # Use bootstrapping to study generate the summary statistics
# # Value is normalized by mean instead of median of controls
# _________

# ## 1 Functions and module

# ### 1.1 Modules

# In[1]:


import pandas as pd
import numpy as np
import copy
import scipy
import random
from scipy.stats import rankdata
import math



#### 2. Bootstrap function
def bootstrap_final_df(focal_df, ref_df, percentile_list, control_gRNA_list, num_replicates, total_gRNA_number, group_variable):
    """
    Performs bootstrapping analysis on focal and reference datasets to generate tumor metrics and relative metrics.

    Parameters:
    - focal_df: DataFrame containing focal group data.
    - ref_df: DataFrame containing reference group data.
    - percentile_list: List of percentiles to calculate tumor metrics.
    - control_gRNA_list: List of control gRNA identifiers.
    - num_replicates: Number of bootstrap replicates to perform.
    - total_gRNA_number: Total number of gRNA to consider in the analysis.
    - group_variable: The column name to group data by (e.g., 'gRNA').

    Returns:
    - DataFrame containing observed and bootstrapped relative tumor metrics.
    """
    # Calculate initial tumor and relative metrics for observed data
    observed_focal_metrics = TumorMetrics(focal_df).calculate_tumor_metrics(percentile_list, [group_variable])
    observed_ref_metrics = TumorMetrics(ref_df).calculate_tumor_metrics(percentile_list, [group_variable])
    observed_relative_metrics = RelativeTumorMetrics(observed_focal_metrics, observed_ref_metrics, group_variable, control_gRNA_list).calculate_relative_tumor_metrics()
    observed_relative_metrics['Bootstrap_id'] = 'Real'

    # Initialize output DataFrame with observed metrics
    output_df = observed_relative_metrics

    # Generate index dictionaries for bootstrapping
    focal_index_dict = generate_index_dictionary(focal_df)
    ref_index_dict = generate_index_dictionary(ref_df)

    # Perform bootstrapping
    for replicate in range(num_replicates):
        focal_bootstrap_indices = nested_bootstrap_index_single(focal_index_dict)
        ref_bootstrap_indices = nested_bootstrap_index_special_single(ref_index_dict, ref_df, total_gRNA_number,group_variable)

        # Generate bootstrap DataFrames
        bootstrap_focal_df = focal_df.loc[focal_bootstrap_indices]
        bootstrap_ref_df = ref_df.loc[ref_bootstrap_indices]

        # Calculate metrics for bootstrap samples
        bootstrap_focal_metrics = TumorMetrics(bootstrap_focal_df).calculate_tumor_metrics(percentile_list, [group_variable])
        bootstrap_ref_metrics = TumorMetrics(bootstrap_ref_df).calculate_tumor_metrics(percentile_list, [group_variable])
        bootstrap_relative_metrics = RelativeTumorMetrics(bootstrap_focal_metrics, bootstrap_ref_metrics, group_variable, control_gRNA_list).calculate_relative_tumor_metrics()
        bootstrap_relative_metrics['Bootstrap_id'] = f'B{replicate}'

        # Append bootstrap metrics to output DataFrame
        output_df = pd.concat([output_df, bootstrap_relative_metrics], ignore_index=True)

    return output_df

def generate_index_dictionary(input_df):
    """Generates a mapping from Sample_ID to DataFrame indices."""
    return {key: group.index.values for key, group in input_df.groupby('Sample_ID')}

def nested_bootstrap_index_single(input_dict):
    """Performs nested bootstrapping, resampling indices within each group."""
    resampled_indices = []
    for group in np.random.choice(list(input_dict.keys()), len(input_dict), replace=True):
        resampled_indices.extend(np.random.choice(input_dict[group], len(input_dict[group]), replace=True))
    return np.array(resampled_indices)

def nested_bootstrap_index_special_single(input_dict, input_df, total_gRNA_number,group_variable):
    """Ensures resampled data contains a minimum number of unique gRNAs."""
    resampled_indices = []
    while len(set(input_df.loc[resampled_indices, group_variable])) < total_gRNA_number:
        resampled_indices = []
        for group in np.random.choice(list(input_dict.keys()), len(input_dict), replace=True):
            resampled_indices.extend(np.random.choice(input_dict[group], len(input_dict[group]), replace=True))
    return np.array(resampled_indices)

#### 3. Tumor Metrics Class
class TumorMetrics:
    def __init__(self, input_df: pd.DataFrame):
        self.raw_df = input_df

    def calculate_tumor_metrics(self, input_percentile: list, group_variable_list: list) -> pd.DataFrame:
        """
        Calculates tumor metrics for each group defined by the group_variable_list.

        :param input_percentile: List of percentiles to calculate.
        :param group_variable_list: List of column names to group the DataFrame by.
        :return: DataFrame with calculated tumor metrics.
        """
        try:
            # Grouping the DataFrame by the specified variables and applying the metric calculation
            temp_df = self.raw_df.groupby(group_variable_list, as_index=False).apply(self._cal_tumor_size_simple, input_percentile)
            return temp_df
        except KeyError as e:
            print(f"Missing column in DataFrame: {e}")
            return pd.DataFrame()

    def _cal_tumor_size_simple(self, x: pd.Series, input_percentile: list) -> pd.Series:
        """
        Auxiliary method to calculate tumor size metrics.

        :param x: Input data for a single group.
        :param input_percentile: List of percentiles to calculate.
        :return: Series with calculated metrics.
        """
        d = {}
        temp_vect = x['Cell_number'].astype(float)
        
        # Measure size
        d['LN_mean'] = self._ln_mean(temp_vect)
        d['Geo_mean'] = self._geometric_mean(temp_vect)
        percentile_values = np.percentile(temp_vect, input_percentile)
        
        for i, percentile in enumerate(input_percentile):
            d[f'{percentile}_percentile'] = percentile_values[i]
        
        d['TTN'] = len(temp_vect)  # Total tumor number
        d['TTB'] = temp_vect.sum()  # Total tumor burden
        return pd.Series(d)

    @staticmethod
    def _ln_mean(input_vector: np.ndarray) -> float:
        """
        Calculates the log-normal mean of the input vector.

        :param input_vector: Input vector of numbers.
        :return: Log-normal mean.
        """
        log_vector = np.log(input_vector)
        mean = log_vector.mean()
        var = log_vector.var() if len(log_vector) > 1 else 0
        return math.exp(mean + 0.5 * var)

    @staticmethod
    def _geometric_mean(input_vector: np.ndarray) -> float:
        """
        Calculates the geometric mean of the input vector.

        :param input_vector: Input vector of numbers.
        :return: Geometric mean.
        """
        log_vector = np.log(input_vector)
        mean = log_vector.mean()
        return math.exp(mean)


#### 4. Relative Tumor Metrics Class
class RelativeTumorMetrics:
    def __init__(self, focal_df: pd.DataFrame, ref_df: pd.DataFrame, group_variable: str, input_control_gRNA_list: list):
        self.focal_df = focal_df
        self.ref_df = ref_df
        self.group_variable = group_variable
        self.control_gRNA_list = input_control_gRNA_list

    def calculate_relative_tumor_metrics(self) -> pd.DataFrame:
        """
        Calculates relative tumor metrics by normalizing metrics in the focal dataset 
        with those in the reference dataset and categorizing each entry.

        :return: DataFrame with calculated relative tumor metrics and annotated sample types.
        """
        # Normalize TTN and TTB
        normalized_metrics_df = self._generate_normalized_metrics(['TTN', 'TTB'],self.group_variable)
        merged_df = self.focal_df.merge(normalized_metrics_df, on=self.group_variable)

        # Calculate relative expression
        self._add_cohort_specific_relative_metrics(merged_df,self.group_variable)

        # Annotate sample type
        merged_df['Type'] = merged_df.apply(
            lambda x: 'Inert' if (x[self.group_variable] in self.control_gRNA_list) else 'Experiment', axis=1
        )
        return merged_df

    def _add_cohort_specific_relative_metrics(self, input_df: pd.DataFrame,group_variable):
        """
        Adds relative metrics for specified columns using the median of control group as the baseline.

        :param input_df: DataFrame with tumor metrics to be normalized.
        """
        control_subset = input_df[input_df[group_variable].isin(self.control_gRNA_list)]
        for column in input_df.columns.difference([group_variable,'Mouse_genotype']):
            relative_column_name = f'{column}_relative'
            input_df[relative_column_name] = input_df[column] / control_subset[column].mean()

    def _generate_normalized_metrics(self, trait_list: list,group_variable) -> pd.DataFrame:
        """
        Generates normalized metrics by dividing metrics in the focal DataFrame by corresponding metrics in the reference DataFrame.

        :param trait_list: List of metrics to be normalized.
        :return: DataFrame with normalized metrics.
        """
        focal_indexed = self.focal_df.set_index(group_variable)
        ref_indexed = self.ref_df.set_index(group_variable).loc[focal_indexed.index]
        normalized_df = pd.DataFrame({group_variable: focal_indexed.index.values})

        for trait in trait_list:
            normalized_column = f'{trait}_normalized'
            normalized_values = np.array(focal_indexed[trait].tolist()) / np.array(ref_indexed[trait].tolist())
            normalized_df[normalized_column] = normalized_values

        return normalized_df


#### 5. Summary Statistics

def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]

def calculate_bootstrap_summary(df, traits_of_interest):
    """
    Calculates summary statistics for bootstrap samples.
    """
    summary = {}
    for trait in traits_of_interest:
        summary.update({
            f'{trait}_95P': df[trait].quantile(0.95),
            f'{trait}_5P': df[trait].quantile(0.05),
            f'{trait}_fraction_greater_than_one': sum(df[trait] > 1) / len(df[trait]),
            f'{trait}_bootstrap_median': df[trait].median(),
            f'{trait}_bootstrap_mean': df[trait].mean(),
            f'{trait}_97.5P': df[trait].quantile(0.975),
            f'{trait}_2.5P': df[trait].quantile(0.025)
        })
    return pd.Series(summary)

def generate_final_summary_dataframe(input_df, traits_of_interest):
    """Generate the final summary DataFrame with FDR-adjusted p-values and two-sided calculations."""
    observed_df = input_df[input_df['Bootstrap_id'] == 'Real']
    bootstrap_df = input_df[input_df['Bootstrap_id'] != 'Real']
    bootstrap_summary = bootstrap_df.groupby('gRNA').apply(calculate_bootstrap_summary, traits_of_interest)
    final_df = observed_df.merge(bootstrap_summary, on='gRNA', suffixes=('', '_bootstrap'))

    for trait in traits_of_interest:
        fraction_col = f'{trait}_fraction_greater_than_one_bootstrap'
        p_value_col = f'{trait}_pvalue'
        fdr_col = f'{trait}_FDR'
        two_side_col = f'{trait}_pvalue_twoside'
        two_side_fdr_col = f'{trait}_FDR_twoside'

        final_df[p_value_col] = final_df[fraction_col].apply(lambda x: min(x, 1-x))
        final_df[two_side_col] = 2 * final_df[p_value_col].clip(upper=1)
        final_df[fdr_col] = fdr(final_df[p_value_col])
        final_df[two_side_fdr_col] = fdr(final_df[two_side_col])

    return final_df

def generate_gene_level_summary_dataframe(input_df, traits_of_interest):
    """Generate gene-level summary DataFrame with FDR-adjusted p-values and two-sided calculations."""
    bootstrap_gene_effect = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(
        ['Targeted_gene_name', 'Bootstrap_id']).apply(calculate_bootstrap_summary, traits_of_interest)
    gene_level_summary = bootstrap_gene_effect.groupby('Targeted_gene_name').mean()
    observed_gene_effect = input_df[input_df['Bootstrap_id'] == 'Real'].groupby('Targeted_gene_name').mean()
    final_gene_level_df = observed_gene_effect.merge(gene_level_summary, on='Targeted_gene_name', suffixes=('', '_bootstrap'))

    for trait in traits_of_interest:
        fraction_col = f'{trait}_fraction_greater_than_one_bootstrap'
        p_value_col = f'{trait}_pvalue'
        fdr_col = f'{trait}_FDR'
        two_side_col = f'{trait}_pvalue_twoside'
        two_side_fdr_col = f'{trait}_FDR_twoside'

        final_gene_level_df[p_value_col] = final_gene_level_df[fraction_col].apply(lambda x: min(x, 1-x))
        final_gene_level_df[two_side_col] = 2 * final_gene_level_df[p_value_col].clip(upper=1)
        final_gene_level_df[fdr_col] = fdr(final_gene_level_df[p_value_col])
        final_gene_level_df[two_side_fdr_col] = fdr(final_gene_level_df[two_side_col])

    return final_gene_level_df

#### 2. Final Summary Statistic for Bootstrap Results
def generate_sgRNA_level_summary_dataframe(input_df, traits_of_interest, group_variable):
    """
    Generates a final summary DataFrame for observed data by integrating bootstrap statistics,
    calculating p-values, FDR adjustments, and their two-sided equivalents for specified traits.

    Parameters:
    - input_df (pd.DataFrame): DataFrame containing both observed data and bootstrap simulations.
    - traits_of_interest (list of str): List of trait names for which the summary statistics are to be calculated.

    Returns:
    - pd.DataFrame: Summary statistics DataFrame for the observed data including bootstrap statistics, p-values, and FDR corrections.
    """
    # Calculate bootstrap statistics for each gRNA
    bootstrap_summary = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(
        group_variable).apply(calculate_bootstrap_summary, traits_of_interest)
    
    # Extract observed data and merge with bootstrap summary statistics
    observed_df = input_df[input_df['Bootstrap_id'] == 'Real'].copy()
    final_summary_df = observed_df.merge(bootstrap_summary, on=group_variable, how='left')

    # Calculate p-values, FDR, and two-sided statistics for each trait of interest
    for trait in traits_of_interest:
        add_statistical_columns(final_summary_df, trait)

    return final_summary_df

def generate_gene_level_summary_dataframe(input_df, traits_of_interest,group_variable):
    """
    Generates a summary DataFrame for gene-level data focusing on specified traits.
    It processes bootstrapped and real data to calculate gene effects and statistical summaries.

    Parameters:
    - input_df: DataFrame with gene-level data, including bootstrap samples.
    - traits_of_interest: List of traits to generate summaries for.

    Returns:
    - DataFrame with combined gene effects for real data and bootstrapping summary statistics.
    """
    # Ensure the DataFrame contains required columns
    required_columns = ['Bootstrap_id', group_variable] + traits_of_interest
    if not all(column in input_df.columns for column in required_columns):
        raise ValueError("Input DataFrame does not contain all required columns.")

    # Process bootstrapped data
    bootstrapped_effects = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(
        [group_variable, 'Bootstrap_id'], as_index=False).apply(
            calculate_combined_gene_effect, traits_of_interest)

    bootstrapping_summary = bootstrapped_effects.groupby(group_variable, as_index=False).apply(
        calculate_bootstrapping_summary, traits_of_interest)

    # Process real data
    real_data_effects = input_df[input_df['Bootstrap_id'] == 'Real'].groupby(
        group_variable, as_index=False).apply(
            calculate_combined_gene_effect, traits_of_interest)

    final_output_df = real_data_effects.merge(bootstrapping_summary, on=group_variable)

    # Calculate additional statistics for each trait
    for trait in traits_of_interest:
        add_statistical_columns(final_output_df, trait)

    return final_output_df

def calculate_combined_gene_effect(group, traits_of_interest):
    """
    Calculates the combined effect of genes for specified traits using weighted averages.

    Parameters:
    - group: Grouped DataFrame segment.
    - traits_of_interest: List of traits to calculate the combined gene effect for.

    Returns:
    - Series with combined effects for specified traits.
    """
    weights = group['TTN_normalized_relative']
    return pd.Series({trait: (group[trait] * weights).sum() / weights.sum() for trait in traits_of_interest})

def calculate_bootstrapping_summary(group, traits_of_interest):
    """
    Calculates summary statistics from bootstrapping results for specified traits.

    Parameters:
    - group: Grouped DataFrame segment from bootstrapped data.
    - traits_of_interest: List of traits to summarize bootstrapping results for.

    Returns:
    - Series with summary statistics for specified traits.
    """
    summary = {}
    for trait in traits_of_interest:
        trait_values = group[trait]
        summary.update({
            f"{trait}_95P": trait_values.quantile(0.95),
            f"{trait}_5P": trait_values.quantile(0.05),
            f"{trait}_fraction_greater_than_one": (trait_values > 1).mean(),
            f"{trait}_bootstrap_median": trait_values.median(),
            f"{trait}_bootstrap_mean": trait_values.mean(),
            f"{trait}_97.5P": trait_values.quantile(0.975),
            f"{trait}_2.5P": trait_values.quantile(0.025)
        })
    return pd.Series(summary)

def add_statistical_columns(df, trait):
    """
    Adds statistical columns to the DataFrame for a given trait.

    Parameters:
    - df: DataFrame to add the statistical columns to.
    - trait: The trait to calculate statistics for.

    Modifies the DataFrame in-place by adding new columns related to p-values and FDR.
    """
    fraction_col = f"{trait}_fraction_greater_than_one"
    pvalue_col = f"{trait}_pvalue"
    fdr_col = f"{pvalue_col}_FDR"
    two_side_col = f"{pvalue_col}_twoside"
    two_side_fdr_col = f"{two_side_col}_FDR"

    df[pvalue_col] = df[fraction_col].apply(lambda x: min(x, 1 - x))
    df[fdr_col] = fdr(df[pvalue_col])
    df[two_side_col] = df[pvalue_col] * 2
    df[two_side_fdr_col] = fdr(df[two_side_col])

def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]
