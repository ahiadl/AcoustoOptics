classdef BackProjection < handle
    %BACKPROJECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vars
        res
        sigMat
    end
    
    methods (Static)
        function uVars = createUserVars()
            uVars.c  = 1480;
            uVars.fs = 40e6;
            uVars.geometry = 'Circular';
            uVars.imageWidth = 30e-3; %[m]
            uVars.mode = 3;
            
            uVars.xSensorPos = []; % horizontal position of the transducer
            uVars.ySensorPos = [];
            uVars.zSensorPos = [];
            uVars.tVec = [];
        end
    end
    
    methods
        function this = BackProjection()
            startupcl;
            this.vars.maxMM = 1e6; %to decide 
        end
        
        function setVars(this, uVars)
            %System Parameters
            this.vars.c    = uVars.c;
            this.vars.fs   = uVars.fs;
            this.vars.tVec = uVars.tVec;
            
            % Time Parameters
            this.vars.timing.offset  = 0; %200e-9; % Light offset from electrical trigger
            this.vars.timing.delay   = 0; % delayed acquisition
%             this.vars.timing.t       = this.vars.timing.offset + ...
%                                        this.vars.timing.delay + ...
%                                        (0:1:2029)*(1/this.vars.fs);
            this.vars.timing.t       = this.vars.tVec;
            
            % Geometry Parameters
            switch uVars.geo
                case 'Circular'
                    r_sensor = (78/2)*10^-3;
                    
                    ang_ini      = -0.7853981633;
                    ang_end      =  3.926990817079;
                    ang_step     =  0.0184799567858;
                    
                    this.vars.geo.angle_sensor = ang_ini:ang_step:ang_end;
                    
                    this.vars.geo.x_sensor = r_sensor*cos(this.vars.geo.angle_sensor); % horizontal position of the transducer
                    this.vars.geo.y_sensor = r_sensor*sin(this.vars.geo.angle_sensor); % vertical position of the transducer
                    this.vars.geo.z_sensor = zeros(length(this.vars.geo.x_sensor),1);

                case 'Linear'
                    
                case 'Custom'
                    this.vars.geo.x_sensor = uVars.xSensorPos; % horizontal position of the transducer
                    this.vars.geo.y_sensor = uVars.ySensorPos;
                    this.vars.geo.z_sensor = uVars.zSensorPos;
                    this.vars.timing.t = uVars.tVec;
                    
                case 'Planar'
                    uVars.xSensorPos = (uVars.xSensorPos-mean(uVars.xSensorPos))*1e-3;
                    uVars.ySensorPos = (uVars.ySensorPos-mean(uVars.ySensorPos))*1e-3;
                    [this.vars.geo.x_sensor,this.vars.geo.y_sensor]=...
                        meshgrid(uVars.xSensorPos,uVars.ySensorPos);
                    this.vars.geo.z_sensor = uVars.zSensorPos*ones(1,length(uVars.xSensorPos)*length(uVars.ySensorPos));
                    this.vars.timing.t = uVars.tVec;
            end
            
            this.vars.geo.sensPosMat  = single([this.vars.geo.x_sensor(:)...
                                                this.vars.geo.y_sensor(:)...
                                                this.vars.geo.z_sensor(:)]);                  
            
            this.vars.numOfSensor = size(this.vars.geo.sensPosMat,1);
            
            % Memory Limitation
            if this.vars.numOfSensor > this.vars.maxMM
                idx = this.vars.maxMM; 
                this.vars.mm = idx;
                while idx < this.vars.numOfSensor
                    idx = idx + this.vars.maxMM;
                    this.vars.mm = [this.vars.mm, idx];
                end
            else
                this.vars.mm = this.vars.numOfSensor;
            end
            
            % Reconstruction Params
            this.vars.recon.dxdy        = this.vars.c/this.vars.fs;
            this.vars.recon.imageWidth  = uVars.imageWidth; 
            this.vars.recon.n           = round(this.vars.recon.imageWidth/this.vars.recon.dxdy);
            this.vars.recon.imageWidth  = this.vars.recon.n * this.vars.recon.dxdy;
            this.vars.recon.mode        = uVars.mode;
            this.vars.recon.reconAxis   = ((-this.vars.recon.n/2+1):1:(this.vars.recon.n/2)) * this.vars.recon.dxdy;  
        end
        
        function vars = getVars(this)
           vars = this.vars; 
        end
        
        function recon = calcBP(this, sigMat)
            this.sigMat = single(permute(squeeze(sigMat),[3,2,1]));
            this.sigMat = reshape(this.sigMat, length(this.vars.tVec), []);
            this.res = this.GPU_BP();
            recon = this.res;
        end
        
        function Recon_cl = GPU_BP(this)
            mm = this.vars.mm;
            Recon_cl = backprojects3d_cl( this.sigMat( :, 1 : mm(1) ),...
                              this.vars.recon.n,...
                              this.vars.geo.sensPosMat( 1 : mm(1), : ),...
                              this.vars.c,...
                              this.vars.recon.mode,...
                              this.vars.timing.t',...
                              this.vars.fs,...
                              this.vars.recon.imageWidth ) ;

            if length(mm)>1
                for i=2:length(mm)
                      pause(0.1);
                      Recon_cl = Recon_cl + backprojects3d_cl( this.sigMat( :, mm(i-1)+1: mm(i) ),...
                                                               this.vars.recon.n,...
                                                               this.vars.geo.sensPosMat( mm(i-1)+1: mm(i), : ),...
                                                               this.vars.c,...
                                                               this.vars.recon.mode,...
                                                               this.vars.timing.t',...
                                                               this.vars.fs,...
                                                               this.vars.recon.imageWidth);
                end
            end
        end
    end
end

