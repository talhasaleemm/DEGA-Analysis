#@title **Define Sample Groups**


# Identify control sample IDs from metadata (water treatment replicates)
control_samples = ["GSM6030870", "GSM6030871", "GSM6030872"]
sample_ids = list(counts_df.columns)
group_labels = ["Control" if sid in control_samples else "SCFA" for sid in sample_ids]
assert len(group_labels) == len(sample_ids), "Mismatch in number of samples and group labels."

# Create a DataFrame for sample metadata
sample_info = pd.DataFrame({'sample': sample_ids, 'group': group_labels})
print(sample_info)
