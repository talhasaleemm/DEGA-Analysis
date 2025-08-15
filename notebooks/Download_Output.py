#@title **Download Output**
import shutil

# Create a ZIP archive of the entire output folder
shutil.make_archive('deg_analysis_output', 'zip', 'deg_analysis_output')
print("Created archive: deg_analysis_output.zip")

# If running in Google Colab, trigger a download dialog
try:
    from google.colab import files
    files.download('deg_analysis_output.zip')
except ImportError:
    print(" To download, locate 'deg_analysis_output.zip' in your file browser and download manually.")

