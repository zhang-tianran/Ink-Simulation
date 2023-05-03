# Script to render a sequence of ply files with Blender. First, setup your scene
# and material for the imported meshes. Scene must be in "OBJECT"-mode.
# Fill in variables in options.
# References.
# See: http://blender.stackexchange.com/questions/24133/modify-obj-after-import-using-python
# and: http://blenderartists.org/forum/showthread.php?320309-How-to-import-ply-files-from-script
import bpy
from typing import Tuple
import numpy as np

# Options.
meshFolder = "/Users/helenhuang/course/cs2240/DaDDi/output"  # Folder without ending "\\".
renderFolder = "/Users/helenhuang/course/cs2240/DaDDi/renders"  # Output folder (without ending "\\").

# meshFolder = "/Users/mandyhe/Documents/Spring2023/Graphics/DaDDi/output"  # Folder without ending "\\".
# renderFolder = "/Users/mandyhe/Documents/Spring2023/Graphics/DaDDi/output"  # Output folder (without ending "\\").

# meshFolder = "output"  # Folder without ending "\\".
# renderFolder = "output"  # Output folder (without ending "\\").
materialName = "Material"  # Material name for the imported object. The Material already needs to be created.
AmountOfNumbers = 1  # Amount of numbers in filepath, e.g., 000010.ply


# Constants.
M_PI = 3.1415926535897932
END_FRAME = 20
ANIM_STEP = 8 # amt of time between frames

# Grid Constants
DIMENSIONS = (14,14,14)
GRID_THICKNESS = 0.01 # thickness of grid lines
BORDER_THICKNESS = 0.1 # thickness of grid border
SHOW_GRID = True # show the grid lines; otherwise, the just the grid border is shown




# define render engine
# bpy.context.scene.render.engine = 'BLENDER_WORKBENCH'
bpy.context.scene.render.engine = 'CYCLES'
#bpy.context.scene.render.engine = 'BLENDER_EEVEE'
#bpy.context.scene.eevee.use_motion_blur = True

# Helper.
def Deg2Rad(degree):
	return degree * (M_PI / 180.0)

def SelectOnlyGivenObject(object):
	# Firs deselect all.
	for iterationObject in bpy.context.scene.objects:
		iterationObject.select_set(False)
	# Then select the given object.
	object.select_set(True)

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
	return folder + "/" + str(frame).zfill(4) + "." + fileEnding

def create_camera(location: Tuple[float, float, float]) -> bpy.types.Object:
    bpy.ops.object.camera_add(location=location)
    return bpy.context.object

def create_light(location):
	# Create light datablock
	light_data = bpy.data.lights.new(name="pointlight-data", type='POINT')
	light_data.energy = 10000
	light_data.distance = 20

	# Create new object, pass the light data 
	light_object = bpy.data.objects.new(name="my-light", object_data=light_data)

	# Link object to collection in context
	bpy.context.collection.objects.link(light_object)

	# Change light position
	light_object.location = location
	return light_object

def create_grid():
	# create the first cube
	bpy.ops.mesh.primitive_cube_add(enter_editmode=False, location=(0,0,0), scale=(0.5,0.5,0.5))
	cube = bpy.context.selected_objects[0]

	for i in range(3):
		mod = cube.modifiers.new('Array_' + str(i), 'ARRAY')
		mod.relative_offset_displace[0] = 0
		mod.relative_offset_displace[i] = 1
		mod.count = DIMENSIONS[i]

	cube.modifiers.new('WD', 'WELD')
	wf = cube.modifiers.new('WF', 'WIREFRAME') 
	wf.thickness = GRID_THICKNESS
	 # shift grid (aka. shift center of first cube after scaling)
	cube.location = [0.5, 0.5, 0.5]

def create_border():
	grid_location = [DIMENSIONS[0]/2, DIMENSIONS[1]/2, DIMENSIONS[2]/2]
	grid_scale = grid_location
	bpy.ops.mesh.primitive_cube_add(enter_editmode=False, location=grid_location, scale=grid_scale)
	cube = bpy.context.selected_objects[0]
	wf = cube.modifiers.new('WF', 'WIREFRAME') 
	wf.thickness = BORDER_THICKNESS


# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------
# CODE FROM: https://github.com/CGArtPython/bpy_building_blocks_examples/tree/main/clean_scene

def purge_orphans():
    if bpy.app.version >= (3, 0, 0):
        bpy.ops.outliner.orphans_purge(
            do_local_ids=True, do_linked_ids=True, do_recursive=True
        )
    else:
        # call purge_orphans() recursively until there are no more orphan data blocks to purge
        result = bpy.ops.outliner.orphans_purge()
        if result.pop() != "CANCELLED":
            purge_orphans()


def clean_scene():
    """
    Removing all of the objects, collection, materials, particles,
    textures, images, curves, meshes, actions, nodes, and worlds from the scene
    """
    if bpy.context.active_object and bpy.context.active_object.mode == "EDIT":
        bpy.ops.object.editmode_toggle()

    for obj in bpy.data.objects:
        obj.hide_set(False)
        obj.hide_select = False
        obj.hide_viewport = False

    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete()

    collection_names = [col.name for col in bpy.data.collections]
    for name in collection_names:
        bpy.data.collections.remove(bpy.data.collections[name])

    # in the case when you modify the world shader
    world_names = [world.name for world in bpy.data.worlds]
    for name in world_names:
        bpy.data.worlds.remove(bpy.data.worlds[name])
    # create a new world data block
    bpy.ops.world.new()
    bpy.context.scene.world = bpy.data.worlds["World"]

    purge_orphans()
    
# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------

# animation code
# referenced: https://blender.stackexchange.com/questions/36902/how-to-keyframe-mesh-vertices-in-python

# def setup_animation(ink_mesh):
# 	action = bpy.data.actions.new("MeshAnimation")
# 	ink_mesh.animation_data_create()
# 	ink_mesh.animation_data.action = action
# 	# instantiate fcurves xyz for each vertex
# 	data_path = "vertices[%d].co"
# 	fcurves = []
# 	for v in ink_mesh.data.vertices:
# 		fcurves.append([action.fcurves.new(data_path % v.index, index =  i) for i in range(3)])
# 	return fcurves

# def insert_keyframe(fcurves, frame, values):
#     for fcu, val in zip(fcurves, values):
#         fcu.keyframe_points.insert(frame, val, options={'FAST'})
	
# def keyframe_vertices(ink_mesh, frame, fcurves):
# 	frame = frame * ANIM_STEP
# 	for i in range(len(ink_mesh.data.vertices)):
# 		curr_vert = ink_mesh.data.vertices[i]
# 		vert_fcurve_xyz = fcurves[i]
# 		insert_keyframe(vert_fcurve_xyz, frame, curr_vert.co)
# 		# fcurves = [action.fcurves.new(data_path % v.index, index =  i) for i in range(3)]
# 		# co_rest = v.co

# 		# for t, value in zip(frames, values):
# 		# 	co_kf = co_rest + value * vec_z
# 		# 	insert_keyframe(fcurves, t, co_kf)






# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------
# ---------- ---------- ---------- ---------- ---------- ----------

def RenderSequence(startFrame = 0, endFrame = 1):
	# Clear the scene:
	clean_scene()

	# create camera
	camera_object = create_camera(location=(0.0, 0.0, 0.0))
	bpy.context.scene.camera = camera_object

	# make the material
	material = bpy.data.materials.new(name="Particle Material")
	material.use_nodes = True
	material_output = material.node_tree.nodes.get('Material Output')
	emission = material.node_tree.nodes.new('ShaderNodeEmission')
	print("new shader node", emission)
	emission.inputs['Strength'].default_value = 15.0
	emission.inputs['Color'].default_value = (1, 0, 0, float(1.0))
	material.node_tree.links.new(material_output.inputs[0], emission.outputs[0])
	light_created = False

	# create grid
	if SHOW_GRID:
		create_grid()
	create_border()


	# enable motion blur
	bpy.context.scene.render.use_motion_blur = True
	# fcurves = None
	ink_mesh = None
	data = []
	


	# Loop over the frames.
	for currentFrame in range(startFrame, endFrame):
		# Import the object (Either obj or ply).
		fullPathToMesh = MeshPath(folder = meshFolder, frame = currentFrame)
		bpy.ops.import_mesh.ply(filepath = fullPathToMesh)

		# Get the just imported object.
		importedObject = bpy.context.object
        
		# Get the smoke material. It has to be named that way. -> did this above in different way
		# material = bpy.data.materials[materialName]

			
		## Camera
		# focus on object
		# for area in bpy.context.screen.areas:
		# 	if area.type == 'VIEW_3D':
		# 		for region in area.regions:
		# 			if region.type == 'WINDOW':
		# 				override = {'area': area, 'region': region}
		# 				bpy.ops.view3d.camera_to_view_selected(override)

		# camera_object.location[2] += 2
		# camera_object.values().lens = 50
		camera_object.location = [7,7,49]

		if ink_mesh is not None:
			ink_coords = np.array([v.co for v in importedObject.data.vertices])
			data.append(ink_coords)
			# coords = [(importedObject.matrix_world @ v.co) for v in importedObject.data.vertices]
			# for i in range(len(ink_mesh.data.vertices)):
			# 	ink_mesh.data.vertices[i].co = coords[i]
			# # ink_mesh.data.vertices = coords
			DeleteObject(importedObject)
		else:
			print("YO")
			# Set the material of the object.
			if len(importedObject.data.materials):
				# assign to 1st material slot
				importedObject.data.materials[0] = material
			else:
				# if there is no material append it
				importedObject.data.materials.append(material)
			## Make and link geometry nodes
			# https://blender.stackexchange.com/questions/259867/geometry-nodes-as-mesh-generation-script
			# 2) Add the GeometryNodes Modifier
			modifier = importedObject.modifiers.new("GeometryNodesNew", "NODES")
			print(modifier.name, modifier)

			# https://blender.stackexchange.com/questions/249763/python-geometry-node-trees/249779#249779
			def new_GeometryNodes_group():
				''' Create a new empty node group that can be used
					in a GeometryNodes modifier.
				'''
				node_group = bpy.data.node_groups.new('GeometryNodes', 'GeometryNodeTree')
				inNode = node_group.nodes.new('NodeGroupInput')
				inNode.outputs.new('NodeSocketGeometry', 'Geometry')
				outNode = node_group.nodes.new('NodeGroupOutput')
				outNode.inputs.new('NodeSocketGeometry', 'Geometry')
				node_group.links.new(inNode.outputs['Geometry'], outNode.inputs['Geometry'])
				# inNode.location = Vector((-1.5*inNode.width, 0))
				# outNode.location = Vector((1.5*outNode.width, 0))
				return node_group
			# In 3.2 Adding the modifier no longer automatically creates a node group.
			# This test could be done with versioning, but this approach is more general
			# in case a later version of Blender goes back to including a node group.
			node_group = None
			if importedObject.modifiers[-1].node_group:
				node_group = importedObject.modifiers[-1].node_group    
			else:
				node_group = new_GeometryNodes_group()
				importedObject.modifiers[-1].node_group = node_group
			nodes = node_group.nodes
			
			meshpoint = nodes.new(type="GeometryNodeMeshToPoints")
			meshpoint.location.x += 400
			meshpoint.location.y -= 50
			# connect
			links = node_group.links
			links.new(nodes["Group Input"].outputs["Geometry"], meshpoint.inputs["Mesh"])
			links.new(meshpoint.outputs["Points"], nodes["Group Output"].inputs["Geometry"])

			# lighting
			if not light_created:
				create_light(camera_object.location)
				light_created = True
			
			ink_mesh = importedObject
			ink_coords = np.array([v.co for v in ink_mesh.data.vertices])
			data.append(ink_coords)

			# setup animation
			# fcurves = setup_animation(ink_mesh)
		
		# bpy.ops.anim.keyframe_insert_menu(type='Location')
		# keyframe_vertices(ink_mesh, currentFrame, fcurves)

		# Render the scene.
		bpy.data.scenes['Scene'].render.filepath = RenderPath(folder = renderFolder, frame = currentFrame)
		bpy.data.scenes["Scene"].camera = camera_object
		bpy.ops.render.render(write_still = True) 

		# Delete the imported object again.
		#DeleteObject(importedObject)
		# DeleteObject(light_object)
	def insert_keyframe(fcurves, frame, values):
		for fcu, val in zip(fcurves, values):
			fcu.keyframe_points.insert(frame, val, options={'FAST'})
	me = ink_mesh.data
	action = bpy.data.actions.new("MeshAnimation")
	me.animation_data_create()
	me.animation_data.action = action

	data_path = "vertices[%d].co"

	frames = list(range(startFrame, endFrame))

	for index, v in enumerate(me.vertices):
		fcurves = [action.fcurves.new(data_path % v.index, index =  i) for i in range(3)]
		for t in frames:
			co_kf = data[t-1][index]
			insert_keyframe(fcurves, t, co_kf)


		
		

# Run the script.
RenderSequence(startFrame = 1, endFrame = END_FRAME)