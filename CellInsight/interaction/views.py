from django.http import HttpResponseRedirect
from django.shortcuts import render

import interaction.forms as forms

def upload_file(request):
    if request.method == 'POST':
        form = forms.UploadFile(request.POST)
        if form.is_valid():
            
            return HttpResponseRedirect()
    else:
        form = forms.UploadFile()
    
    return render(request, 'upload.html', {'form':form})
