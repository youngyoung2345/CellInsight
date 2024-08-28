from django import forms

class UploadFileForm(forms.Form):
    file = forms.FileField()
    file_format = forms.ChoiceField(choices=[('h5ad', 'H5AD'),('h5', 'H5'),('txt', 'Text/Clustering File'),('csv', 'CSV File')])

class QCForm(forms.Form):
    min_counts = forms.IntegerField(label='Min Counts', initial=2000)
    min_genes = forms.IntegerField(label='Min Genes', initial=500)
    max_genes = forms.IntegerField(label='Max Genes', initial=7000)
    pct_counts_mt = forms.FloatField(label='Pct Counts MT', initial=20.0)

class GeneSearchForm(forms.Form):
    gene_name = forms.CharField(label='Gene Name', max_length=100, required=True)