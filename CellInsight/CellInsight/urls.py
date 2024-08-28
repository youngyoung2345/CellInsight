"""
URL configuration for CellInsight project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/5.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.contrib import admin
from django.urls import path
from interaction import views
from django.conf import settings
from django.conf.urls.static import static
urlpatterns = [
    path('admin/', admin.site.urls),
   # path('upload/', views.upload_file, name='upload_file'),
   # path('success/', views.success, name='success'),

    path('', views.only_render,{'html': 'welcome.html'},name='welcome'),
    path('preprocessing/', views.preprocessing, name='preprocessing'),
    path('qc_process/', views.qc_process, name='qc_process'),
    path('mapcell_process/', views.mapcell_process, name='mapcell_process'),
    path('umap/', views.umap_view, name='umap_view'),
    path('markersearch/', views.markersearch, name='markersearch'), 
    path('genesearch/', views.genesearch, name='genesearch'),

   path('search/', views.only_render, {'html': 'search.html'}, name='search'),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
