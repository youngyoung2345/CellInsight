from django import forms

class UploadFileForm(forms.Form):
    file = forms.FileField()
    file_format = forms.ChoiceField(choices=[('h5ad', 'H5AD'),('txt', 'Text/Clustering File'),('cellranger', 'CellRanger Output'),('csv', 'CSV File')])


class QCForm(forms.Form):
    min_counts = forms.IntegerField(label='Min Counts', initial=2000)
    min_genes = forms.IntegerField(label='Min Genes', initial=500)
    max_genes = forms.IntegerField(label='Max Genes', initial=7000)
    pct_counts_mt = forms.FloatField(label='Pct Counts MT', initial=20.0)