import sys
# # kratos_libs="/home/enric/kratos/libs"
# # kratos_applications="/home/enric/kratos/applications"
# # kratos_scripts="/home/enric/kratos/kratos/python_scripts"
# # sys.path.append(kratos_libs)
# # sys.path.append(kratos_applications) 
# # sys.path.append(kratos_scripts) 
import os.path
kratos_libs=os.path.join(os.path.dirname(__file__),'../libs')
kratos_applications=os.path.join(os.path.dirname(__file__),'../applications')
kratos_scripts=os.path.join(os.path.dirname(__file__),'../kratos/python_scripts')
sys.path.append(kratos_libs)
sys.path.append(kratos_applications) 
sys.path.append(kratos_scripts) 
