function u0 = TahoeIC(x) % Initial condition for GHF
global Xrel0 Xrelf IntTime iteration tspan sol
idx = find(tspan <= x, 1, 'last' );
if iteration == 1
    u0 = Xrel0 + (Xrelf-Xrel0)*sin(pi*x/(2*IntTime));
else
    u0 = squeeze(sol(end,idx,:));
end
end
% --------------------------------------------------------------------------