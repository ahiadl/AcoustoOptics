classdef stages < handle
    %STAGES supports up to 3 stages labeled ['X' 'Y' 'Z']
    
    properties
        hw
        type
        
        connected
        hardwareAvailable
        
        vars  
        curPos % according to log Idx
        
        axOrder            % arrange according to user request with relation to log ID
        logId    %1   2   3 % logical id arranged according to logical Order (fixed)
        hwId     %4   5   6 % hw id arranged according to hw order.
        
        % Logical id is always 1:numOfStages.
        % hw id is the order in which the hw loaded the stages. it is set
        % when connecting to hw.
        % The only thing the user can affect is how each axes is called
        % (one letter) so later operation can be mafe by either logical id
        % or ax name. hw id is transparent to user as it may be different
        % from one hw to other.
        
    end
    
    methods
        function this = stages(type)
            this.type = type;
            
            switch type
                case 'Zaber'
                    this.hw = zaberStages('COM3');
                case 'PI'
                    this.hw = PIStages();
                otherwise
                    this.hw = [];
                    fprintf("STAGES: Can't find the requested hw.\n")
                    return;
            end
            
            this.vars.axOrderName = ['X', 'Y', 'Z'];
            this.connected = false;
            this.hardwareAvailable = false;
        end

        function connect(this)
            this.connected = this.hw.connect();

            if this.connected 
                this.vars.numOfStages = this.hw.getNumOfStages();
                this.vars.validAxes   = zeros(1, this.vars.numOfStages);
                this.curPos = zeros(1,this.vars.numOfStages);
                
                % Fill logic stage number with hw stage number:
                % e.g. logId(1) = 1 -> hwId(1) = 2; 
                this.logId = 1:this.vars.numOfStages;
                [axDef, hwIdxDef] = this.hw.getDefaultAssign();
                this.hwId = hwIdxDef;
                
                this.assignStagesAxes(axDef);
                                
                this.connected = true;
                this.hardwareAvailable = true;
            else
                this.connected = false;
                this.hardwareAvailable = false;
                fprintf("STAGES: Can't connect to stages\n")
            end
        end
        
        % Set Functions
        function assignStagesAxes(this, axOrder) 
            % Function aligns the axName vector and the id Vector
            % so user may refer to axes as 'X' or 'Y' instead of 1, 2.
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't assign\n");
               return;
            end
            
            axIn = length(axOrder);
            ax = unique(axOrder);
            
            if axIn > length(ax) || axIn ~= this.vars.numOfStages
                fprintf("Stages: Can't assign stages. Please check the requested assignment.\n")
                return
            end

            this.axOrder = '';
            for i = 1:axIn
                this.axOrder(i) = axOrder(i);
            end
            
            this.updatePosition();
            this.vars.fullRange = this.getFullRange(); 
            [ax, up, low] = this.hw.getDefaultLimits();
            this.setLimits(ax, up, low);
        end

        function setLimits(this, ax, upper, lower)
            for i=1:length(ax)
                idx = this.getLogIdFromAxName(ax(i));
                this.vars.limits.upper(idx) = upper(i);
                this.vars.limits.lower(idx) = lower(i);
            end
        end
        
        function setVelocity(this, vel)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, set velocity\n");
               return;
            end
            
        end
        
        % Get Functions
        function range = getFullRange(this)
            % This function returns each stage full possible range.
            % sorted by axOrder
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't get range\n");
               return;
            end
            
            range = zeros(1,this.vars.numOfStages);
            for ax = 1:this.vars.numOfStages
                curId = this.hwId(ax);
                range(ax) = this.hw.readPosition(curId);
            end
        end

        function updatePosition(this)
            % This function updates the curPos field by directly reading
            % each stage position. result is ordered by axOrder
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't update position\n");
               return;
            end
            
            this.vars.curPos = [];
            for i = 1:this.vars.numOfStages
                curId = this.hwId(i);
                this.curPos(i) = this.hw.readPosition(curId);
            end
        end
        
        function curLogId = getLogIdFromAxName(this, ax)
            % this function translate ax name (e.g. 'X') to its idx 
            % in axOrder (e.g. 1)
            curLogId = zeros(1,length(ax));
            for j = 1:length(ax)
                for i =1:length(this.axOrder)
                    if strcmp(this.axOrder(i), ax(j))
                        curLogId(j) = i; 
                    end
                end
            end
        end

        function hwId = getHwIdFromLogId(this, id)
            % this function translate stage id (e.g. 2) to its hwId 
            hwId = zeros(1,length(id));
            for i = 1:length(id)
                hwId(i) = this.hwId(id(i));
            end
        end
        
        function id = getLogicIDs(this)
            % Sorted by axes
            id = this.logId;
        end
        
        function pos = getPosVar(this)
            pos = this.curPos;
        end
        
        function pos = getPosition(this)
            this.updatePosition();
            pos = this.getPosVar;
        end
        
        % Move Functions
        function moveRelId(this, id, rel)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't move\n");
               return;
            end
            
            idPos = this.curPos(id);
            if idPos + rel > this.vars.limits.upper(id)
                fprintf("STAGES: can't move stage %s - upper Limit Exception\n", this.axOrder(id))
                return;
            end
            
            if idPos + rel < this.vars.limits.lower(id)
                fprintf("STAGES: can't move stage %s - lower Limit Exception\n", this.axOrder(id))
                return;
            end
            
            curHwId = this.getHwIdFromLogId(id);
            this.hw.moveRel(curHwId, rel)
            this.updatePosition();
        end
        
        function moveAbsId(this, id, pos)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't move\n");
               return;
            end
            
            if pos > this.vars.limits.upper(id)
                fprintf("STAGES: can't move stage %s - upper Limit Exception\n", this.axOrder(id))
                return;
            end
            if pos < this.vars.limits.lower(id)
                fprintf("STAGES: can't move stage %s - lower Limit Exception\n", this.axOrder(id))
                return;
            end
            
            curHwId = this.getHwIdFromLogId(id);
            this.hw.moveAbs(curHwId, pos);
            this.updatePosition();
        end

        function moveRelAx(this, ax, rel)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't move\n");
               return;
            end
            
            curLogId = this.getLogIdFromAxName(ax);
            this.moveRelId(curLogId, rel)
        end
        
        function moveAbsAx(this, ax, pos)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't move\n");
               return;
            end
            
            curLogId = this.getLogIdFromAxName(ax);
            this.moveAbsId(curLogId, pos);
        end
        
        %non blocking
        function moveNBRelId(this, id, rel)
        
        end
        
        function moveNBAbsId(this, id, rel)
        
        end
        
        function moveNBRelAx(this, id, rel)
        
        end
        
        function moveNBAbsAx(this, id, rel)
        
        end
        
        function home(this, idax)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't go home\n");
               return;
            end
            
            if sum(this.vars.limits.lower > 0) % TODO: implement lower bound that is not zero
                fprintf("STAGES: Home is out of limits\n");
                return;
            end
            if ischar(idax)
                curLogId = this.getLogIdFromAxName(idax);
                curHwId = this.getHwIdFromLogId(curLogId);
                this.hw.home(curHwId)
            else
                curHwId = this.getHwIdFromLogId(idax);
                this.hw.home(curHwId);
            end
        end
        
        function allHome(this)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't move home\n");
               return;
            end
            
            if sum(this.vars.limits.lower > 0) % TODO: implement lower bound that is not zero
                printf("STAGES: Home is out of limits\n");
                return;
            end
            this.hw.homeAll();
        end

        function moveStagesForIdentification(this, id)
            if ~(this.hardwareAvailable && this.connected)
               fprintf("STAGES: Stages not available, can't identify by moving\n");
               return;
            end
            
            span = 10;
            
            curHwId = this.getHwIdFromLogId(id);
            
            if (this.curPos(curHwId) + span) > this.vars.limits.upper(curHwId)
                this.hw.moveRelId(id, -span);
                this.hw.moveRelId(id,  span);
            elseif this.curPos(curHwId) - span < this.vars.limits.lower(curHwId)
                this.hw.moveRelId(id,   span);
                this.hw.moveRelId(id, - span);
            else
                fprintf("STAGES: span is too big for assignment test\n");
            end
        end
    end
end

