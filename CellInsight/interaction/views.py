from django.http import HttpResponseRedirect
from django.shortcuts import render, redirect
from django.utils.http import urlencode
from preprocessing.prepo import preprocess_and_visualize
import interaction.forms as forms
import search.marker_search as marker_search
import pandas as pd
import os
import preprocessing.Umap as Umap
def handle_uploaded_file(f):
    media_path = 'media'
    if not os.path.exists(media_path):
        os.makedirs(media_path)
    file_path = os.path.join(media_path, f.name)
    with open(file_path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    return file_path

def preprocessing(request):
    if request.method == 'POST':
        form = forms.UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            file = request.FILES['file']
            file_format = form.cleaned_data['file_format']
            file_path = handle_uploaded_file(file)
            qc_form = forms.QCForm()
            return render(request, 'preprocessing.html', {
                'file_uploaded': True,
                'form': form,
                'qc_form': qc_form,
                'file_path': file_path,
                'file_format': file_format
            })
    else:
        form = forms.UploadFileForm()
        qc_form = forms.QCForm()
        return render(request, 'preprocessing.html', {
            'file_uploaded': False,
            'form': form,
            'qc_form': qc_form
        })


def qc_process(request):
    if request.method == 'POST':
        qc_form = forms.QCForm(request.POST)
        file_path = request.POST.get('file_path')
        file_format = request.POST.get('file_format')

        if qc_form.is_valid():
            min_counts = qc_form.cleaned_data['min_counts']
            min_genes = qc_form.cleaned_data['min_genes']
            max_genes = qc_form.cleaned_data['max_genes']
            pct_counts_mt = qc_form.cleaned_data['pct_counts_mt']

            # Violin Plot만 생성
            violin_plot, _ = preprocess_and_visualize(file_path, file_format, min_counts, min_genes, max_genes, pct_counts_mt)

            return render(request, 'preprocessing.html', {
                'violin_plot': violin_plot,
                'qc_form': qc_form,
                'show_mapcell_button': True,
                'file_uploaded': True,
                'file_path': file_path,
                'file_format': file_format,
                'min_counts': min_counts,
                'min_genes': min_genes,
                'max_genes': max_genes,
                'pct_counts_mt': pct_counts_mt
            })

    return HttpResponseRedirect('/preprocessing/')

def welcome(request):
    return render(request, 'welcome.html')

def mapcell_process(request):
    if request.method == 'POST':
        file_path = request.POST.get('file_path')
        file_format = request.POST.get('file_format')
        min_counts = float(request.POST.get('min_counts'))
        min_genes = float(request.POST.get('min_genes'))
        max_genes = float(request.POST.get('max_genes'))
        pct_counts_mt = float(request.POST.get('pct_counts_mt'))

        # Normalization 및 UMAP Plot 생성
        _, umap_plot = preprocess_and_visualize(file_path, file_format, min_counts, min_genes, max_genes, pct_counts_mt)

        # UMAP 플롯 페이지로 리디렉션
        umap_plot = umap_plot.replace('\\', '/')
        query_params = urlencode({'umap_plot': umap_plot})

        # GET 파라미터로 전달하여 UMAP 플롯 페이지로 리디렉션
        return redirect(f'/umap/?{query_params}')

    return HttpResponseRedirect('/preprocessing/')

def umap_view(request):
    user_umap_plot = request.GET.get('umap_plot')
    s3_umap_plot = None

    if request.method == 'POST':
        folder_name = request.POST.get('s3_file')
        delimiter = request.POST.get('delimiter', '\t')

        if folder_name:
            try:
                # 클러스터 파일 가져오기 및 UMAP 생성
                cluster_file2 = Umap.fetch_cluster_files(folder_name)
                
                if cluster_file2:
                    # 선택된 클러스터 파일로 UMAP 생성
                    s3_umap_plot = Umap.fetch_and_process_file(cluster_file2, delimiter)
                else:
                    print("No cluster files found.")
            except Exception as e:
                print(f"Error fetching and processing file {folder_name}: {str(e)}")
                return render(request, 'umap.html', {
                    'user_umap_plot': user_umap_plot,
                    's3_umap_plot': None,
                    's3_files': Umap.fetch_s3_folder_list(),
                    'error_message': f"Error processing file: {str(e)}"
                })

    # S3 파일 목록 가져오기
    s3_files = Umap.fetch_s3_folder_list()

    return render(request, 'umap.html', {
        'user_umap_plot': user_umap_plot,
        's3_umap_plot': s3_umap_plot,
        's3_files': s3_files,
    })

    
def markersearch(request):
    # marker_search.py에서 DataFrame 데이터를 가져옴
    markergene_df = marker_search.fetch_markers()

    # organ별로 cell type 리스트 생성
    organ_cell_types = markergene_df.groupby('organ')['cell type'].apply(list).to_dict()
    for organ in organ_cell_types:
        cell_type_counts = pd.Series(organ_cell_types[organ]).value_counts()
        organ_cell_types[organ] = [f"{cell_type} ({count})" for cell_type, count in cell_type_counts.items()]

    # items()의 결과를 리스트로 변환하여 템플릿에 전달
    organ_cell_types_items = list(organ_cell_types.items())

    # DataFrame을 HTML 테이블로 변환
    markergene_html = markergene_df.to_html(classes='table table-striped', index=False)

    return render(request, 'markersearch.html', {
        'organ_cell_types_items': organ_cell_types_items,
        'markergene_html': markergene_html,
    })

def search(request):
    return render(request, 'search.html')
