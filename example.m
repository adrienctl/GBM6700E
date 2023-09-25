%% Usage examples for the function GetObject()

load Calib_Beads2D.mat
load Calib_Beads3D.mat
load Vertebrae2D.mat

global PA0 PA20;  % values to define the two views
PA0 = 1;
PA20 = 2;

% Cell array containing beads from both views
Beads2D = cell(2,1);
Beads2D{1} = Beads2D_PA0;
Beads2D{2} = Beads2D_PA20;

%% Retrieve 2D coords of a bead available in both views

bead_2D = zeros(2,2);
name = 'A_1_1';
for viewpt = PA0:PA20,
    obj = GetObject(Beads2D{viewpt},name);
    bead_2D(viewpt,:) = obj.coord;
end
clear obj

%% Retrieve 2D coords in the 2 views of the points forming a vertebra

name = 'Vertebra_T1';
vertT1_PA0  = GetObject(Vertebrae_PA0,name);
vertT1_PA20 = GetObject(Vertebrae_PA20,name);
vertebra_PA0  = zeros(6,2);
vertebra_PA20 = zeros(6,2);
for k = 1:6,
    vertebra_PA0(k,:)  = vertT1_PA0.points2D(k).coord;
    vertebra_PA20(k,:) = vertT1_PA20.points2D(k).coord;
end
clear k

%% Retrieve the 3D coords of a set of beads by their names

beads = ['A_1_1';'A_2_4';'A_3_4';'A_4_3';'A_5_5';...
          'B_1_1';'B_1_5';'B_3_3';'B_4_2';'B_5_3'];
beads_3D = zeros(length(beads),3);
for k = 1:length(beads)
    obj = GetObject(Calib_Beads3D, beads(k,:));
    beads_3D(k,:) =  obj.coord;
end
clear obj k

