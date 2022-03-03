classdef DAQengineNew < handle
%% Handheld Optoacoustic Tomography acquisition and reconstruction engine.
%  Use this class together with the .NET4 DAQ-dlls and OpenCL recon dlls.
%  This class provides easy access to data acquisition, reconstruction,
%  post-processing and multi-spectral analysis.


properties ( Access = protected )
    
    % Fetcher Thread Sync Param
    ftsp = struct( 'locked', 0, ...
        'sigMat', [], ...
        'sigMat_buffer', [], ...
        'channels', 0, ...
        'samples', 0, ... 
        'samples_delay', 0, ...
        'connected_channels', [], ...
        'num_av', 0, ...
        'NumberOfFrames', 1, ...
        'preview', 0,...
        'recon_counter', 0,...
        'recon_rate', 0,...
        'preview_mode', 0,...
        'frame_counter', 1, ...
        'av_counter', 0, ...
        'live_acq', 0, ...
        'live_maxFrames', 3000, ...
        'thread', [], ...
        'daq_success', 0, ...
        'wbh', 0, ...
        'imh', [], ...
        'init_success', 0) ;
    acuiring;
    
end

%
%
%
properties ( GetAccess = public )
    
    % reconstruction parameters
    conf = [] ;
    
end

%
%
%
methods ( Access = public )
        
    function HE = COLDengine()
        
        HE.conf = [];
        
    end
        
    %
    %   initiate DAQ
    %
    function success = initDAQ( this, samples, channels, mac, connected_channels, samples_delay )
        
        if( isempty( connected_channels ) )
            
            this.ftsp.connected_channels = 1:channels ;
            
        else
            
            this.ftsp.connected_channels = connected_channels ;
            
        end
        
        this.ftsp.samples           = samples ;
        this.ftsp.samples_delay     = samples_delay ;
        this.ftsp.channels          = channels ;
        disp(['MaxFrames: ', num2str(this.ftsp.live_maxFrames)]);
        this.ftsp.sigMat            = [] ;
        this.ftsp.sigMat            = zeros( samples, length( this.ftsp.connected_channels ), this.ftsp.live_maxFrames ) ;
        this.ftsp.sigMat_buffer     = zeros( samples, length( this.ftsp.connected_channels ) ) ;

        NET.addAssembly( [ cd '\Peripherals\DAQFiles\fmDAQWrapper.dll'] ) ;
        tw = fmDAQWrapper.ThreadWrapper( channels, mac ) ;
        this.ftsp.thread = tw ;
        
        % locked check
        disp(['locked = ', num2str(this.ftsp.locked)]);
        addlistener( tw, 'DAQEvent', @this.daemon_callback ) ;
        this.ftsp.wbh = waitbar( 0, 'Initializing DAQ. Please wait ...', 'name', 'DAQ Initialization' ) ;
        tw.Main() ;
        waitfor( this.ftsp.wbh ) ;
        success = this.ftsp.init_success ;
        %this.ftsp.wbh = 0 ;
        
    end

    %
    %   free DAQ
    %
    function success = closeDAQ( this )
        
        try
            
            this.ftsp.thread.stopThread() ;
            delete( this.ftsp.thread ) ;
            fprintf( '\nDAQ closed.\n' ) ;
            success = 1 ;
            
        catch
            
            success = 0 ;
            
        end
        
    end
        
    %
    %   lock DAQ and initiate data accquisition
    %
    function lockDAQ( this, lock ) %imh )
        
        this.resetData();
        this.ftsp.locked = lock ;
                        
    end
    
    %
    %   swtich live preview on/off
    %
    function LivePreview( this, preview, conf, rate, imh, mode )

        this.conf = [];
        this.conf = conf;
        this.ftsp.recon_rate = rate;
        this.ftsp.recon_counter = 1;
        this.ftsp.imh = imh ;
        this.ftsp.preview_mode = mode;
        this.ftsp.preview = preview;
           
    end
    
    %
    %
    %
    function sigMat = DoAcquisition( this, num_av, num_frames )
        
        if num_frames <= this.ftsp.live_maxFrames
            
            this.ftsp.NumberOfFrames = num_frames ;
            this.ftsp.num_av         = num_av ; 

            %   reset internal variables
            this.resetData();
            this.ftsp.frame_counter  = 1 ;
            this.ftsp.av_counter     = 1 ;
            this.ftsp.daq_success    = 0 ;
            this.ftsp.live_acq       = 1;
            
            this.acuiring = true;
%             this.ftsp.wbh            = waitbar( 0, 'Acquisition running, please wait ...' ) ;

           while(this.acuiring)
               pause(0.001);
           end
            %this.ftsp.wbh = 0 ;
            %this.ftsp.live_acq = 0;

            sigMat                  = this.ftsp.sigMat(:, :, :) ;
            
        else
           
            disp(['Max number of frames exceeded: ', num2str(this.ftsp.live_maxFrames)]);
            
        end
            
    end
    
    %
    %
    %
    function RunningAcquisition( this, num_av, num_frames, hh )
        
        if num_frames <= this.ftsp.live_maxFrames
            
            this.ftsp.NumberOfFrames = num_frames + 3;
            this.ftsp.num_av         = num_av ; 

            %   reset internal variables
            this.resetData();
            this.ftsp.frame_counter  = 1 ;
            this.ftsp.av_counter     = 1 ;
            this.ftsp.daq_success    = 0 ;
            this.ftsp.live_acq       = 1;
            
            this.ftsp.wbh            = hh; %waitbar( 0, 'Acquisition running, please wait ...' ) ;
       
        else
           
            disp(['Max number of frames exceeded: ', num2str(this.ftsp.live_maxFrames)]);
            
        end
            
    end 
    
    function sigMat = GetData( this )
        
        sigMat = this.ftsp.sigMat(:, :, :) ;
        
    end
         
    %
    %   reset sigMat
    %
    function resetData( this )

        this.ftsp.sigMat = [];
        this.ftsp.sigMat = zeros( this.ftsp.samples, length( this.ftsp.connected_channels ), this.ftsp.NumberOfFrames );
        
        this.ftsp.sigMat_buffer = [];
        this.ftsp.sigMat_buffer = zeros( this.ftsp.samples, length( this.ftsp.connected_channels ) ) ;

    end
                                                  
end
    
methods ( Access = private )
    
    %
    %   fetcher thread
    %
    function daemon_callback( this, source, args )

        %%   DAQ error routine
        if( args.fmError )

            disp( ['Error!'] ) ;
            fprintf( ['\n' char( args.pleoraMsg ) '\n'] ) ;
            fprintf( ['\n' char( args.errorMsg ) '\n'] ) ;
            this.ftsp.daq_success = 0 ;

        end

        %%   DAQ timeout routine
        if( args.fmTimeout )

            disp( ['Timeout!'] ) ; 

        end

        %%   standard routine
        if( this.ftsp.locked && ~args.fmError && ~args.fmTimeout )
            % locked check
            disp(['locked = ', num2str(this.ftsp.locked)]);
            disp( ['Event!' ]); 
            
            tmp_aux = args.Data.double;
%             size(tmp_aux)
            tmp_aux2 = tmp_aux(1:(this.ftsp.samples+2)*this.ftsp.channels);
%             size(tmp_aux2)
            tmp = reshape( tmp_aux2, this.ftsp.samples+2, this.ftsp.channels ) ;
%             tmp = reshape( args.Data.double, this.ftsp.samples, this.ftsp.channels ) ;

            if( this.ftsp.samples_delay )

                tmp = [ zeros( this.ftsp.samples_delay, length( this.ftsp.connected_channels ) ) ;
                tmp( 3:this.ftsp.samples+2, this.ftsp.connected_channels ) ] ;
                tmp( this.ftsp.samples+1:end, : ) = [] ;

            else

                tmp = tmp( 3:this.ftsp.samples+2, this.ftsp.connected_channels ) ;

            end
            
            %%  Live Preview
            if ( this.ftsp.preview )

                if( ~isempty( this.ftsp.imh ) ) && this.ftsp.recon_counter == this.ftsp.recon_rate

                    this.ftsp.recon_counter = 1;

                    if  this.ftsp.preview_mode == 1

                        [temp, Recon] = Eat(this.conf, -tmp);
                        set( this.ftsp.imh(1, 1),'CData', squeeze(max(Recon(:,:,:), [], 1)) ) ;
                        set( this.ftsp.imh(1, 2),'CData', squeeze(max(Recon(:,:,:), [], 2)) ) ;
                        set( this.ftsp.imh(1, 3),'CData', squeeze(max(Recon(:,:,:), [], 3)) ) ;

                    elseif this.ftsp.preview_mode == 2

                        [temp, Recon] = Eat(this.conf, tmp);
                        set( this.ftsp.imh(1, 1),'CData', squeeze(max(Recon(:,:,:), [], 1)) ) ;
                        set( this.ftsp.imh(1, 2),'CData', squeeze(max(Recon(:,:,:), [], 2)) ) ;
                        set( this.ftsp.imh(1, 3),'CData', squeeze(max(Recon(:,:,:), [], 3)) ) ;

                    elseif this.ftsp.preview_mode == 3

                        imagesc(tmp);
                        title('Raw Signals');
                        xlabel('Detectors');
                        ylabel('Time Steps');
                        colorbar();

                    elseif this.ftsp.preview_mode == 4

                        plot(tmp(:, 1));
                        title(['Raw Signal of Detector 1']);
                        xlabel('Time Step');
                        ylabel('a.u.');

                    end

                else
                    
                    this.ftsp.recon_counter = this.ftsp.recon_counter + 1;
                    
                end
                
            end
                        
            %%  DAQ acquistion management
            if( this.ftsp.live_acq && tmp(1, 1) ~= 0 )
                
                 %%  Do what ever you want with tmp ... inset code here
           
                 % ... bla bla bla
                 %   DoRecon(tmp, imh);
                 % [temp, Recon] = Eat(this.conf, sigMat);
                 
                
                %   acquistion running
                if this.ftsp.av_counter <= this.ftsp.num_av
                    
                    this.ftsp.sigMat_buffer = this.ftsp.sigMat_buffer + tmp ;
                    disp(['Frame: ', num2str(this.ftsp.frame_counter), '/', num2str(this.ftsp.NumberOfFrames), ', Average: ', num2str(this.ftsp.av_counter), '/', num2str(this.ftsp.num_av)]);
                    this.ftsp.av_counter = this.ftsp.av_counter + 1;

                else

                    this.ftsp.sigMat( :, :, this.ftsp.frame_counter ) = this.ftsp.sigMat_buffer ./ this.ftsp.num_av ;
                    this.ftsp.sigMat_buffer = [];
                    this.ftsp.sigMat_buffer = zeros( this.ftsp.samples, length( this.ftsp.connected_channels ) ) ;

                    
                    %   save frame frame_counter + 1   
                    this.ftsp.sigMat_buffer = this.ftsp.sigMat_buffer + tmp ;
                    this.ftsp.frame_counter = this.ftsp.frame_counter + 1 ;
                    this.ftsp.av_counter = 2;

                    %   stop acquisition when number of frames are reached
                    if this.ftsp.frame_counter > this.ftsp.NumberOfFrames
                    
                        this.ftsp.live_acq = 0;  
                        this.ftsp.locked = 0;
                        this.ftsp.daq_success = 1;
                        
                    else
                        
                        disp(['Frame: ', num2str(this.ftsp.frame_counter), '/', num2str(this.ftsp.NumberOfFrames), ', Average: ', num2str(this.ftsp.av_counter-1), '/', num2str(this.ftsp.num_av)]);
  
                        
                    end
                     
                end                             
                   
            end
            
        end
        if( ~this.ftsp.locked )
                
            this.ftsp.init_success = ~args.fmError ;
%             delete( this.ftsp.wbh ) ;
            this.acuiring = false;
        end 
    end

end

end

