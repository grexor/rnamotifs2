import rnamotifs2
import os
import glob

def init():
    rnamotifs2.path.root_folder = os.path.abspath(os.path.join(os.path.abspath(__file__), "..", ".."))
    rnamotifs2.path.comps_folder = os.path.join(root_folder, "comps")
