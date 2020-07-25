function h = plotcube(L,origin,Alpha_index,color,DCM)
vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];

l = L(1); w = L(2); h = L(3);
% Height change
vert(5,3) = h+vert(5,3); vert(6,3) = h+vert(6,3); vert(7,3) = h+vert(7,3); vert(8,3) = h+vert(8,3); 

% With change
vert(3,2) = w+vert(3,2); vert(4,2) = w+vert(4,2); vert(8,2) = w+vert(8,2); vert(7,2) = w+vert(7,2);

% lenght chage
vert(2,1) = l+vert(2,1); vert(3,1) = l+vert(3,1); vert(6,1) = l+vert(6,1); vert(7,1) = l+vert(7,1);

% Rotation
vert_transpose = DCM*vert'; vert = vert_transpose';

hold on
translation_matrix = repmat(origin,size(vert,1),1);
vert = vert + translation_matrix;

h = patch('Vertices',vert,'Faces',fac,...
      'FaceVertexCData',color,'FaceColor','flat');
alpha(Alpha_index) 
% view(3)
% axis vis3d

end
