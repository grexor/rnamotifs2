import rnamotifs2
import os

# comps_folder is the directory that holds the per-comparison subdirectories
# (comps_folder/<name>/<name>.tab, <name>.config, r1s/, rnamap/, ...).
#
# Resolution order:
#   1. an explicit rnamotifs2.path.set_comps_folder(...) call
#   2. the RNAMOTIFS2_COMPS environment variable
#   3. <current working directory>/comps  (what run_example.sh relies on)
#
# expressRNA (or any embedding application) should call set_comps_folder()
# once after importing rnamotifs2.


def init():
    if getattr(rnamotifs2.path, "comps_folder", None) is not None:
        return
    env = os.environ.get("RNAMOTIFS2_COMPS")
    if env:
        rnamotifs2.path.comps_folder = os.path.abspath(os.path.expanduser(env))
    else:
        rnamotifs2.path.comps_folder = os.path.join(os.getcwd(), "comps")
    # kept for backwards compatibility; nothing in the package uses it
    rnamotifs2.path.root_folder = os.path.dirname(rnamotifs2.path.comps_folder)


def set_comps_folder(folder):
    rnamotifs2.path.comps_folder = os.path.abspath(os.path.expanduser(folder))
    rnamotifs2.path.root_folder = os.path.dirname(rnamotifs2.path.comps_folder)
