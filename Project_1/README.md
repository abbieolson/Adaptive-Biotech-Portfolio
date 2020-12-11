## Project 1 - Classifier identifying and analyzing potential inaccuracies

#### *Languages and Tools:*
* Python
  * pandas
  * NumPy
  * SciPy
  * scikit-learn
  * seaborn
  * Matplotlib
  * Plotly
  * Jupyter
* HTML
----------
### Pre-Processing:
* pandas
* seaborn
```python3
# Make Dataframe and Plot Missing Values
df = pd.read_csv('df.csv', low_memory=False)
df = pd.DataFrame(df, index=df.index, columns=df.columns)

df = df[(df.thing1 == "other_thing")]
df = df[(df.thing2 == "another_thing")]

target = 'target'

# plot missing values
ans = df.drop(target, axis=1).isnull().sum().sort_values(ascending=False)
plt.figure(figsize=(12,2))
sns.heatmap(pd.DataFrame(data=ans[ans>0], columns=['Missing Values']).T, annot=True, cbar=False, cmap='viridis', annot_kws={'rotation': 90})

plt.savefig('missing_vals.png')
```

### Data Visualizations:
* seaborn
* Matplotlib
* NumPy
```python3
# Correlation Heatmap
def make_pair_plot(df):
    '''Function to make correlation plots with seaborn'''
    sns.set(style="white")

    mask = np.triu(np.ones_like(df_corr, dtype=np.bool))

    mask[np.triu_indices_from(mask)] = True

    # matplotlib figure setup
    f, ax = plt.subplots(figsize=(12, 12))

    # custom colormap
    cmap = sns.diverging_palette(250, 10, as_cmap=True)

    # draw heatmap with aspect ratio
    sns.heatmap(df_corr, mask=mask, cmap=cmap, vmax=.3, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5})

    # save with matplotlib
    plt.savefig('corr_plot.png')
```
* Plotly
```python3
# Interactive Bar Chart
fig5 = make_subplots(rows=2, cols=1)

fig5.append_trace(go.Bar(
    x=list(v_amp_yes_gene_in_ct['GeneName']),
    y=list(v_amp_yes_gene_in_ct['Count']),
    name='"All" (Amp)',
    marker_color='rgb(95, 182, 239)',
    width=0.4
), row=1, col=1)

fig5.append_trace(go.Bar(
    x=list(v_amp_no_gene_in_ct['GeneName']),
    y=list(v_amp_no_gene_in_ct['Count']),
    name='"All" (No Amp)',
    marker_color='rgb(204, 204, 255)',
    width=0.4
), row=1, col=1)

fig5.append_trace(go.Bar(
    x=list(v_amp_yes_gene_ct['GeneName']),
    y=list(v_amp_yes_gene_ct['Count']),
    name='"In" (Amp)',
    marker_color='rgb(168, 241, 115)',
    width=0.4
), row=2, col=1)

fig5.append_trace(go.Bar(
    x=list(v_amp_no_gene_ct['GeneName']),
    y=list(v_amp_no_gene_ct['Count']),
    name='"In" (No Amp)',
    marker_color='rgb(255, 176, 140)',
), row=2, col=1)

fig5.update_layout(
    # template='plotly_dark',
    legend_title_text='Status',
    title="Fig. 5: Gene Counts",
    yaxis_title="Count",
    font=dict(
        size=14,
    ))

# change the bar mode
fig5.update_layout(barmode='group')

pio.write_image(fig5, 'fig5.jpeg', scale=20)
fig5.show()
```

### Random Forest
* scikit-learn
* Matplotlib
```python3
# Selecting important features, creating a Random Forest model, and assessing with a ROC curve
df = df
target = 'target'

numeric_features = df.select_dtypes(include=['int64', 'float64']).drop(target, axis=1).columns
categorical_features = df.select_dtypes(include=['object']).columns

X = df.drop([target, 'GC', 'AT', 'nt_length', 'homopol_count'], axis=1)
y = df[target]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, stratify=df[target], random_state = 100, shuffle=True)

# exhaustive grid search
model_RF = RandomForestClassifier()
n_estimators = [10, 50, 100]
max_features = ['sqrt', 'log2']
# define grid search
grid = dict(n_estimators=n_estimators,max_features=max_features)
cv = RepeatedStratifiedKFold(n_splits=10, n_repeats=3, random_state=1)
grid_search = GridSearchCV(estimator=model_RF, param_grid=grid, cv=cv, scoring='roc_auc')
grid_result_RF = grid_search.fit(X, y)
# summarize results
print("Best from ALL Features: %f using %s" % (grid_result_RF.best_score_, grid_result_RF.best_params_))
means = grid_result_RF.cv_results_['mean_test_score']
stds = grid_result_RF.cv_results_['std_test_score']
params = grid_result_RF.cv_results_['params']
for mean, stdev, param in zip(means, stds, params):
    print("%f (%f) with: %r" % (mean, stdev, param))

# best model from grid search
RF = RandomForestClassifier(**grid_result_RF.best_params_)
RF.fit(X_train, y_train)

y_pred_RF=RF.predict(X_test)
y_pred_RF_train=RF.predict(X_train)

# get top 20 features
feat_importances = pd.Series(RF.feature_importances_, index=X.columns)
feat_importances.nlargest(19).plot(kind='barh')
plt.savefig('imp_feats.png')

# subset to top 15 features
feature_list = list(X)
indices = np.argsort(feat_importances)[::-1]
ranked_index_2 = [feature_list[i] for i in indices]

imp_feats = ranked_index_2[0:16]

# AUROC
rfc_pred = rfc.predict_proba(X_test)
false_positive_rateRF, true_positive_rateRF, thresholdRF = roc_curve(y_test, rfc_pred[:,1])
roc_aucRF = auc(false_positive_rateRF, true_positive_rateRF)
plt.figure(figsize = (8,8))
plt.title('Receiver Operating Characteristic')
plt.plot(false_positive_rateRF, true_positive_rateRF, color = 'green', label = 'RF AUC = %0.2f' % roc_aucRF)
plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1], linestyle = '--')
plt.axis('tight')
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.savefig('rf_roc_auc.png')
```
### Write Report
* HTML
* Python
```html
HTML script for writing final report

html_string = '''
<html>
    <head>
      <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
      <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
      <style>body{ margin:0 100; background:Gainsboro; }</style>

      <style>
      img {
      display: block;
      margin-left: auto;
      margin-right: auto;}
      </style>

      </head>
      <body>
      <p style="text-align:center;font-size:38px;">Redacted</p>

      <p style="text-align:center;font-size:25px;">Data</p>

      <p style="text-align:center;font-size:20px;">Raw Data</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      <div style="overflow-x:auto;">
      ''' + summary_table_1 + '''
      </div>

      <p style="text-align:center;font-size:20px;">Transformed Data</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      <div style="overflow-x:auto;">
      ''' + summary_table_2 + '''
      </div>


      <p style="text-align:center;font-size:25px;">Exploratory Analysis</p>

      <p style="text-align:center;font-size:25px;">Figure 1</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig1_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 2</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig2_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 3</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig3_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 4</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig4_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 5</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig5_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 6</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig6_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

      <p style="text-align:center;font-size:25px;">Figure 7</p>
      <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
      ''' + fig7_html + '''
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>

    <p style="text-align:center;font-size:25px;">Build a Random Forest Classifier</p>

        <p style="text-align:center;font-size:25px;">Correlation Plot of Data</p>
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
    <img src="../PROJECTS/corr_plot.png" style="width:50%">
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>


    <p style="text-align:center;font-size:25px;">Top Twenty Features from Classifier</p>
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
    <img src="../PROJECTS/rf_imp_feats.png" style="width:50%">
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>


    <p style="text-align:center;font-size:25px;">AUROC of Random Forest Classifier</p>
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>
    <img src="../PROJECTS/rf_roc_auc.png" style="width:50%">
    <p style="text-align:center;font-size:16px;">Stuff about the figure.</p>


    </body>
</html>'''

```
```python3
# write selected results to HTML filename
with open("proj_report.html", 'w') as f:
    f.write(html_string)
```
