classdef DiffuseSimulator
    %DIFFUSESIMULATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function this = DiffuseSimulator()
            toastThreadCount(0);
        end
        
        
        function createMeshName(numOfLayers, meshWidth, layersDepths)
            
             meshName = sprintf("%dLayerSlabMesh-%dx%d-DxW-mm",  vars.meshDepth, vars.meshWidth);
            
            
        end
        
        function outputArg = createMesh(obj,meshVars)
            meshPath
            
            
            switch vars.geometry
                case 'rectangle'
                    meshName = this.createMeashName()
                   
                    meshPath = "./Meshes";
                    meshFilename = sprintf("%s.msh", meshName);
                    meshFullName = sprintf("%s/%s", meshPath, meshFilename);
                    resDirName = sprintf("E:/AcoustoOptics/Simulations/Results/%dx%d-DxW-mm", meshDepth(j), meshWidth(k));
                    
                    %Create the mesh from the geometrical points defined in the .geo file
                    fprintf("Constructing Mesh\n"); 
                    if ~exist(meshFullName, 'file')
                        fprintf("Can't locate the mesh. Generating...\n"); 
                        str(1) = sprintf("depth = %d;", vars.meshDepth);
                        str(2) = sprintf("width = %d;", vars.meshWidth);
                        replaceLineInFile("./1LayerBox.geo", [1,2],  str);
                        system (sprintf("gmsh -3 ./1LayerBox.geo -o %s", meshFullName));
                    else
                        fprintf("Located the mesh. Loading without generating.\n");
                    end
                    
                    % Load the mesh into toastMesh object.
                    % Notice: If the .msh file was already constructed the above code is
                    % redundant and the mesh should only be loaded with the code below.
                    fprintf("Loading Mesh\n");
                    tic;
                    mesh = toastMesh(meshFullName,'gmsh');
                    toc

                    if ~exist(resDirName, 'dir')
                        mkdir(resDirName)
                    end
                    
                case 'ball'
                    
                case 'other'
            end
        end
        
        
        
    end
end

