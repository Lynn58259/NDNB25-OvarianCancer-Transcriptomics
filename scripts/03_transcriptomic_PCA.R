python<-Sys.getenv('PYTHON',unset=Sys.which('python'))
if(!nzchar(python))stop('Python with NumPy and pandas is required; set PYTHON to its executable.')
status<-system2(python,shQuote(normalizePath('scripts/transcriptomic_pca.py')))
if(status!=0)stop('Original NumPy PCA computation failed')
