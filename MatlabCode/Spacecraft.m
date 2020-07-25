classdef Spacecraft    
% -------------------------------------------------------------------------
%  This class descripbe the physical properties of a solar sail
%
% [inputs]: - [Array] : 7 by 1 vector composed of [L l cb M m rhos rhod Cr]
%
% [outputs]: - [I_total] : Total moment of inertia of the sailcraft
% -------------------------------------------------------------------------

   properties
     L    % Length of the solar sail
     l    % width of the solar sail     
     cb   % Control boom
     M    % Mass of the bus
     m    % Mass of the coner supports
     rhos % rate of specular reflection
     rhod % rate of diffusion
     Cr   % reflectivity coefficient
     AMratio % Area/mass ratio
     It % Transversal moment of inertial
     In % Axial moment of inertia
   end
   methods
       
       function Current_Spacraft = Spacecraft(Array)
           if length(Array) == 9
               Current_Spacraft.L         = Array(1);
               Current_Spacraft.l         = Array(2);
               Current_Spacraft.cb        = Array(3);
               Current_Spacraft.M         = Array(4);
               Current_Spacraft.m         = Array(5);
               Current_Spacraft.rhos      = Array(6);
               Current_Spacraft.rhod      = Array(7);
               Current_Spacraft.Cr        = 1+Array(6);
               Current_Spacraft.AMratio   = Array(1)*Array(2)/(Array(4)+4*Array(5));
               Current_Spacraft.It        = Array(8);
               Current_Spacraft.In        = Array(9);
           end
           
       end
       
      function I_total = Moment_Of_Inertia_Calculator(obj)
          
         I_total = [obj.It obj.It obj.In]';
         
      end
      
   end
end















