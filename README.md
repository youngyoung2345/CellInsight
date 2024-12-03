# CellInsight

Single cell analysis 효율 증진 툴 개발 프로젝트

<br>

## 프로젝트 목적

Single Cell RNA-Seq(scRNA-seq, 단일 세포 시퀀싱)이란?
: 생명체의 유전체를 단일 세포 수준에서 시퀀싱하는 생명공학 기술  
*시퀀싱 : DNA의 염기서열(A, T, G, C)을 읽는 것

**Annotation**이란? 세포 데이터의 각 클러스터가 어떤 세포 종류인지 labeling하는 작업  
❗자동화 도구가 있더라도 결국 사람이 레퍼런스 데이터를 확인하며 최종 판단 필요  

![annotation](https://github.com/user-attachments/assets/8a7ece3d-9aa1-4e09-a4d9-d658b6b7b6f6)

✔️Annotation 작업의 효율성을 향상시키고 보다 편리한 분석 환경을 제공하기 위한 웹 플랫폼 구현

<br>

## 연구 방법

**1. Single Cell 데이터베이스 통합**  
현재 가장 많이 사용되고 있는 PanglaoDB와 Single Cell Portal의 2가지 데이터베이스를 통합하여 제공한다. 데이터를 통합할 뿐만 아니라 각 데이터베이스가 제공하는 조금씩 다른 기능들도 통합하여 제공함으로써 여러 사이트를 이용할 필요가 없도록 한다.  

**2. 레퍼런스 데이터 추천 기능 구현**  
사용자가 원하는 조건의 레퍼런스 데이터를 효율적으로 찾기 위해 컨텐츠 기반 필터링 알고리즘을 이용하여 추천시스템을 구현한다.  

**3. 효율적인 분석 환경 제공**  
Single Cell Analysis를 직접 수행해본 사용자 경험을 기반으로 하여 가장 효율적인 화면 구성 및 기능 배치를 제공한다.  

<br>

## 연구 과정

![architecture](https://github.com/user-attachments/assets/06e5d8d6-24e8-4c9f-a0ca-5cfdd9ddbeea)

**1. 데이터 다운로드 및 전처리**  
- PanglaoDB  
(1) bulk data 다운로드 -> RData 형식의 파일에서 sparse matrix 추출  
(2) scRNA-seq 샘플 데이터 크롤링 -> sparse matrix의 각 feature에 맞게 매핑

- Single Cell Portal  
(1) curl을 이용하여 스터디별 데이터 수동 다운로드  
(2) mtx, h5/h5ad, csv 등 파일 형식별로 읽어들여 데이터 조립  
-> 형식별로 데이터 조립 방식 상이, 스터디마다 파일 형식이 다름  

**2. 통합 데이터베이스 및 웹 플랫폼 구축**  
- 네이버 클라우드 플랫폼 버킷에 전처리 끝난 데이터 적재  
- 레퍼런스 데이터를 불러와서 시각화 결과를 볼 수 있는 웹 플랫폼 제작  
-> Django 이용, 버킷 연동  

**3. 레퍼런스 데이터 검색 및 추천 시스템 구현**  
- Search 기능 구현 : 스터디 검색, Gene Marker 검색, Gene 검색  
- 레퍼런스 데이터 추천시스템 구현 : PanglaoDB, Single Cell Portal의 서로 다른 데이터 형식에 맞게 TF-IDF를 이용한 컨텐츠 기반 필터링 적용  

<br>

## 연구 결론
![overview](https://github.com/user-attachments/assets/0759e2cd-a66b-4a2a-abbf-0fc85bf97014)
![home](https://github.com/user-attachments/assets/aed733a7-5726-482c-b242-8db6c62fc7e0)
![search](https://github.com/user-attachments/assets/b3bc2b6f-7deb-41bf-aae4-2807534d426e)
![recommendation](https://github.com/user-attachments/assets/81305f88-0acc-4e24-a596-5206d0adc350)

CellInsight 웹사이트 시연 영상  
https://carpal-money-c50.notion.site/CellInsight-81bee13926e84ff392a0f4deeac23bc9?pvs=4

<br>

## Members  

| 류여진 | 이하영 |
| :-: | :-: |
| <img src='https://avatars.githubusercontent.com/u/88676496?v=4' height=130 width=130></img> | <img src='https://avatars.githubusercontent.com/u/134286859?v=4' height=130 width=130></img> |
| <a href="https://github.com/ryj8075" target="_blank"><img src="https://img.shields.io/badge/GitHub-black.svg?&style=round&logo=github"/></a> | <a href="https://github.com/youngyoung2345" target="_blank"><img src="https://img.shields.io/badge/GitHub-black.svg?&style=round&logo=github"/></a> |
