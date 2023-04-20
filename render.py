# Script to render a sequence of ply files with Blender. First, setup your scene
# and material for the imported meshes. Scene must be in "OBJECT"-mode.
# Fill in variables in options.
# References.
# See: http://blender.stackexchange.com/questions/24133/modify-obj-after-import-using-python
# and: http://blenderartists.org/forum/showthread.php?320309-How-to-import-ply-files-from-script
import bpy
from typing import Tuple

# Options.
meshFolder = "output"  # Folder without ending "\\".
renderFolder = "output"  # Output folder (without ending "\\").
materialName = "Material"  # Material name for the imported object. The Material already needs to be created.
AmountOfNumbers = 1  # Amount of numbers in filepath, e.g., 000010.ply

# Constants.
M_PI = 3.1415926535897932

# Helper.
def Deg2Rad(degree):
	return degree * (M_PI / 180.0)

def SelectOnlyGivenObject(object):
	# Firs deselect all.
	for iterationObject in bpy.context.scene.objects:
		iterationObject.select = False
	# Then select the given object.
	object.select = True

# Delete an object, see: http://blender.stackexchange.com/questions/27234/python-how-to-completely-remove-an-object
def DeleteObject(object):
	# Cache the currrent mode we are in.
	#oldMode = bpy.context.mode
	# Set it to object mode.
	#bpy.ops.object.mode_set(mode = 'OBJECT')
	# Select only the given object.
	SelectOnlyGivenObject(object)
	# Delete the object and set the mode back to where it was.
	bpy.ops.object.delete()
	#bpy.ops.object.mode_set(mode = oldMode)

def MeshPath(folder = "", frame = 0, fileEnding = "ply"):
	return folder + "/" + str(frame).zfill(AmountOfNumbers) + "." + fileEnding

def RenderPath(folder = "", frame = 0, fileEnding = "png"):
	return folder + "/" + str(frame).zfill(AmountOfNumbers) + "." + fileEnding

def create_camera(location: Tuple[float, float, float]) -> bpy.types.Object:
    bpy.ops.object.camera_add(location=location)
    return bpy.context.object

def set_camera_params(camera: bpy.types.Camera,
                      focus_target_object: bpy.types.Object,
                      lens: float = 85.0,
                      fstop: float = 1.4) -> None:
    # Simulate Sony's FE 85mm F1.4 GM
    camera.sensor_fit = 'HORIZONTAL'
    camera.sensor_width = 36.0
    camera.sensor_height = 24.0
    camera.lens = lens
    camera.dof.use_dof = True
    camera.dof.focus_object = focus_target_object
    camera.dof.aperture_fstop = fstop
    camera.dof.aperture_blades = 11
    
def add_track_to_constraint(camera_object: bpy.types.Object, track_to_target_object: bpy.types.Object) -> None:
    constraint = camera_object.constraints.new(type='TRACK_TO')
    constraint.target = track_to_target_object
    constraint.track_axis = 'TRACK_NEGATIVE_Z'
    constraint.up_axis = 'UP_Y'

def RenderSequence(startFrame = 0, endFrame = 1):
	# Loop over the frames.
	for currentFrame in range(startFrame, endFrame):
		# Import the object (Either obj or ply).
		fullPathToMesh = MeshPath(folder = meshFolder, frame = currentFrame)
		# bpy.ops.import_scene.obj(filepath = full_path_to_file)
		bpy.ops.import_mesh.ply(filepath = fullPathToMesh)

		# Get the just imported object.
		importedObject = bpy.context.object

		# Set its orientation. We need to do this,
		# as PreonLab meshes have another up-axis.
		importedObject.rotation_euler = (Deg2Rad(90), Deg2Rad(0), Deg2Rad(0))
        
		# Get the smoke material. It has to be named that way.
		material = bpy.data.materials[materialName]
        # Get material
		# Set the material of the object.
		if len(importedObject.data.materials):
			# assign to 1st material slot
			importedObject.data.materials[0] = material
		else:
			# if there is no material append it
			importedObject.data.materials.append(material)
			
		## Camera
		camera_object = create_camera(location=(0.0, 0.0, 0.0))

		add_track_to_constraint(camera_object, importedObject)
		set_camera_params(camera_object.data, importedObject, lens=72, fstop=0.5)

		# Render the scene.
		bpy.data.scenes['Scene'].render.filepath = RenderPath(folder = renderFolder, frame = currentFrame)
		bpy.data.scenes["Scene"].camera = camera_object
		bpy.ops.render.render(write_still = True) 

		# Delete the imported object again.
		#DeleteObject(importedObject)

# Run the script.
RenderSequence(startFrame = 1, endFrame = 2)