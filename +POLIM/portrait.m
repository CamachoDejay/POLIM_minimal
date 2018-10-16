classdef portrait < handle
    %PORTRAIT polarization portrait, once it knows its intensity and angles
    %it can calculate modulation depths and also display data
    
    properties (SetAccess = private)
        I_ex_em
        I_ex
        I_em
        exAngRad
        emAngRad
        Mex
        Pex
        Mem
        Pem
        LS
    end
    
    methods
        function obj = portrait(I, exAngRad, emAngRad)
            %PORTRAIT Construct an instance of this class
            %   Detailed explanation goes here
            
            assert(ismatrix(I), 'Intensity musb be a matrix')
            assert(size(I,1) == length(exAngRad), 'Ex angles dont match intensity matrix size');
            assert(size(I,2) == length(emAngRad), 'Em angles dont match intensity matrix size');
            
            obj.I_ex_em = I;
            obj.exAngRad = exAngRad;
            obj.emAngRad = emAngRad;
            
        end
        
        function [exAngV, emAngV, intV] = linearize(obj)
            
            % Vectorizing the angles
            [EmGrid, ExGrid] = meshgrid(obj.emAngRad, obj.exAngRad);
            % I want angles to be a one row vector size = [1xn].
            exAngV = ExGrid(:)';
            emAngV = EmGrid(:)';

            % Vectorizing the portrait
            portLength = length(exAngV);
            intV = reshape(obj.I_ex_em,portLength ,1);
            
        end
        
        function set.exAngRad(obj,rowVec)
            assert(isrow(rowVec), 'ex angles must be a row vector')
            assert(min(rowVec) >= 0, 'ex angles must be >= 0')
            assert(max(rowVec) < pi, 'ex angles must be in radians 0 <= ex < pi')
            
            obj.exAngRad = rowVec;
           
        end
        
        function set.emAngRad(obj,rowVec)
            assert(isrow(rowVec), 'em angles must be a row vector')
            assert(min(rowVec) >= 0, 'em angles must be >= 0')
            assert(max(rowVec) < pi, 'em angles must be in radians 0 <= ex < pi')
            
            obj.emAngRad = rowVec;
           
        end
        
        function getModulations(obj)
            obj.getMexPex;
            obj.getMemPem;
            
            % Luminescence shift
            LumShift = obj.Pex - obj.Pem; 
            if LumShift>pi/2
                LumShift=LumShift-pi;
            elseif LumShift<-pi/2
                LumShift=LumShift+pi;
            end
            obj.LS = LumShift;

        end
        
        function getMexPex(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % integration over emission angles
            obj.I_ex = sum(obj.I_ex_em,2);
            % normalization
            I_exNorm = obj.I_ex ./ sum(obj.I_ex);
            % fit
            [Y, ~] = POLIM.fitModulation(obj.exAngRad',I_exNorm,pi/180);
            % parcing output
            obj.Mex = Y(1);
            obj.Pex = Y(2);
        end
        
        function getMemPem(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % integration over emission angles
            obj.I_em = sum(obj.I_ex_em,1)';
            % normalization
            I_emNorm = obj.I_em ./ sum(obj.I_em);
            % fit
            [Y, ~] = POLIM.fitModulation(obj.emAngRad',I_emNorm,pi/180);
            % parcing output
            obj.Mem = Y(1);
            obj.Pem = Y(2);
        end
        
        function showPortrait(obj)
            
            figure()
            contourf(obj.I_ex_em)
            axis image
            a = gca;
            
            tmp = a.XTick;
            tmp = cellstr(num2str(tmp'));
            a.XTickLabel = tmp;
            
            tmp = a.YTick;
            tmp = cellstr(num2str(tmp'));
            a.YTickLabel = tmp;
            colormap 'gray'
            
            xlabel('Excitation Polarization Angle')
            ylabel('Emission Polarization Angle')
            a.FontSize = 14;
            
            colorbar
        end
    end
end

