classdef DAQ < handle
    %DAQ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        core
        vars
    end
    
    methods
        function this = DAQ()
            this.vars.mac = '00-11-1c-02-04-00';
            this.vars.NumberOfChannels = 256;
            this.vars.fs = 40e6;
            
            this.vars.uVars.numOfSamples = 0;
            this.vars.uVars.activeChannels = zeros(1,256);
            this.vars.uVars.delay = 0;
            
            this.vars.changedParams = true;
            this.core = DAQengine();
            this.vars.connected = false;
        end
        
        function configAndConnect(this)
            if this.vars.changedParams || ~(this.vars.connected)
                if this.vars.connected
                    this.disconnect();
                    this.vars.connected = false;
                end

                this.core.initDAQ( this.vars.numOfSamples,...
                                   this.vars.NumberOfChannels, ...
                                   this.vars.mac,...
                                   this.vars.activeChannels,...
                                   this.vars.delay );
                this.vars.connected = true;
            end
        end
        
        function setUserVars(this, uVars)
            % HW vars
            this.vars.prevuVars = this.vars.uVars;
            this.vars.uVars     = uVars;
            
%             this.vars.inputImpedance = uVars.inputImpedance;
            this.vars.activeChannels  = uVars.activeChannels;
            
            % Timing Vars
            this.vars.delay           = uVars.delay;
            this.vars.numOfSamples    = uVars.numOfSamples;
            
            % Acquisition Vars
            this.vars.numOfAvg    = uVars.numOfAvg;
            this.vars.numOfFrames = uVars.numOfFrames;
            
            this.vars.numOfActiveChannels = sum(this.vars.activeChannels>0);
            %TODO: % set impedance, frequency.
            
            diff1 = (this.vars.uVars.numOfSamples ~= this.vars.prevuVars.numOfSamples);
            diff2 = (sum(this.vars.uVars.activeChannels ~= this.vars.prevuVars.activeChannels) > 0);
            diff3 = (this.vars.uVars.delay ~= this.vars.prevuVars.delay);
            
            this.vars.changedParams  = (diff1 + diff2 + diff3) > 0;
        end
        
        function vars = getVars(this)
            vars = this.vars;
        end

        function res = acquire(this)
            this.core.lockDAQ(1);
            res = this.core.DoAcquisition(this.vars.numOfAvg, this.vars.numOfFrames);
            this.core.lockDAQ(0);
        end
        
        function disconnect(this)
            this.core.lockDAQ(0);
            this.core.closeDAQ();
            this.vars.connected = false;
        end
        
        function resetDAQ(this)
            if this.vars.connected
                this.disconnect();
            end
            this.core = DAQengine();
            this.configAndConnect()
        end
    end
end

