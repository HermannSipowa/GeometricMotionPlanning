function u0 = TahoeIC(x) % Initial condition for GHF
global Xrel0 Xrelf IntTime
u0 = Xrel0 + (Xrelf-Xrel0)*sin(pi*x/(2*IntTime));
end
% --------------------------------------------------------------------------