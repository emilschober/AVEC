import os
import runpy

# Ensure relative paths resolve from this directory
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE_DIR)

# Execute the main application file without triggering its __main__ block
globals_dict = runpy.run_path("20251110_Tool_With_Changes_Marlen.py", run_name="avec_app")

app = globals_dict.get("app")
if app is None:
    raise RuntimeError("WSGI: Flask 'app' not found after loading 20251110_Tool_With_Changes_Marlen.py")

