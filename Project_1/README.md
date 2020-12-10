## Project 1 - Classifier identifying and analyzing potential causes of inaccuracy
Languages and Tools:
* Python
  * pandas
  * NumPy
  * SciPy
  * scikit-learn
  * seaborn
  * Plotly
* HTML

Data Visualizations:

```python3
# seaborn

def make_pair_plot(df_corr):
    '''Function for correlation plots'''
    sns.set(style="white")

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(df_corr, dtype=np.bool))
    mask[np.triu_indices_from(mask)] = True

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(12, 12))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(250, 10, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(df_corr, mask=mask, cmap=cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})
    plt.savefig('corr_plot.png')
```
