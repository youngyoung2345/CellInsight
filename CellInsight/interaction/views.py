from django.http import HttpResponseRedirect, JsonResponse
from django.shortcuts import render, redirect
from django.utils.http import urlencode

import interaction.forms as forms
import search.marker_search as marker_search
import pandas as pd
import os
import preprocessing.umap as Umap
import search.marker_search as marker_search
import preprocessing.PanglaoDB_proc_python as PanglaoDB_proc_python
from preprocessing.proc import * 
def only_render(request, html):
    return render(request, html)

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

def handle_uploaded_file(f):
    media_path = 'media'
    if not os.path.exists(media_path):
        os.makedirs(media_path)
    file_path = os.path.join(media_path, f.name)
    with open(file_path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
    return file_path


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

            data = load_data(file_path, file_format)
            preprocessed_data = preprocess_data(data, min_counts, min_genes, max_genes, pct_counts_mt)
            violin_plot = draw_and_save_violin_plot(preprocessed_data)

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

def mapcell_process(request):
    if request.method == 'POST':
        file_path = request.POST.get('file_path')
        file_format = request.POST.get('file_format')
        min_counts = float(request.POST.get('min_counts'))
        min_genes = float(request.POST.get('min_genes'))
        max_genes = float(request.POST.get('max_genes'))
        pct_counts_mt = float(request.POST.get('pct_counts_mt'))

        # Normalization 및 UMAP Plot 생성

        data = load_data(file_path, file_format)
        preprocessed_data = preprocess_data(data, min_counts, min_genes, max_genes, pct_counts_mt)
        normalized_data = normalize_data(preprocessed_data)
        umap_plot = draw_and_save_umap_plot(normalized_data)

        # UMAP 플롯 페이지로 리디렉션
        umap_plot = umap_plot.replace('\\', '/')
        query_params = urlencode({'umap_plot': umap_plot})

        # GET 파라미터로 전달하여 UMAP 플롯 페이지로 리디렉션
        return redirect(f'/umap/?{query_params}')

    return HttpResponseRedirect('/preprocessing/')

def umap_view(request):
    user_umap_plot = request.GET.get('umap_plot')
    s3_umap_plot = None
    panglao_umap_plot = None
    csv_data = None  # 초기화
    if request.method == 'POST':
        folder_name = request.POST.get('s3_file')
        delimiter = request.POST.get('delimiter', '\t')

        if folder_name:
            try:
                if 'PanglaoDB' in folder_name:
                    # PanglaoDB 파일을 선택한 경우 처리
                    panglao_umap_plot, obs_data, uns_data = Umap.fetch_and_process_file_Panglao(folder_name)
                    csv_data = PanglaoDB_proc_python.import_additional_data(folder_name)
                else:
                    # Singlecellportal 파일을 선택한 경우 처리
                    cluster_file2 = Umap.fetch_cluster_files(folder_name)
                    
                    if cluster_file2:
                        s3_umap_plot = Umap.fetch_and_process_file(cluster_file2, delimiter)
                    else:
                        print("No cluster files found.")
            except Exception as e:
                print(f"Error fetching and processing file {folder_name}: {str(e)}")
                return render(request, 'umap.html', {
                    'user_umap_plot': user_umap_plot,
                    's3_umap_plot': None,
                    'panglao_umap_plot': None,
                    's3_files': Umap.fetch_s3_folder_list(),
                    'csv_data': csv_data,
                    'error_message': f"Error processing file: {str(e)}"
                })

    # S3 파일 목록 가져오기
    s3_files = Umap.fetch_s3_folder_list()
    panglao_files = Umap.fetch_s3_folder_list_Panglao()

    # singlecellportal 및 panglaoDB 파일들을 하나의 리스트로 합침
    combined_files = s3_files + [(folder, folder) for folder in panglao_files]
    csv_data_list = csv_data.to_dict(orient='records') if csv_data is not None else []
    return render(request, 'umap.html', {
        'user_umap_plot': user_umap_plot,
        's3_umap_plot': s3_umap_plot,
        'panglao_umap_plot': panglao_umap_plot,
        's3_files': combined_files,  # 통합된 파일 리스트를 전달
        'csv_data': csv_data_list,
        'obs_data': obs_data if 'obs_data' in locals() else None,
        'uns_data': uns_data if 'uns_data' in locals() else None,
    })

    
from django.shortcuts import render
import search.marker_search as marker_search

def markersearch(request):
    selected_organ = request.GET.get('organ', 'All')
    selected_cell_type = request.GET.get('cell_type', 'All')

    # 필터링된 결과를 가져옴
    organ_cell_types, html_code = marker_search.fetch_markers(selected_organ, selected_cell_type)

    # 선택된 organ에 해당하는 cell types 리스트 가져오기
    cell_types = organ_cell_types.get(selected_organ, []) if selected_organ != 'All' else []

    # 템플릿으로 데이터 전달
    return render(request, 'markersearch.html', {
        'organ_cell_types_items': organ_cell_types.items(),
        'selected_organ': selected_organ,
        'selected_cell_type': selected_cell_type,
        'cell_types': cell_types,
        'html_code': html_code,
    })

def search(request):
    return render(request, 'search.html')