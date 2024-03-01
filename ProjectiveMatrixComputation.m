
function ProjectiveMatrixComputation

% Importing the images into mathlab
f = imread("b1.jpg");
g = imread("b2.jpg");

%Selecting points on the images
figure(1), imshow(f);
[Image1X,Image1Y] = ginput(8);
hold on;
plot(Image1X,Image1Y); hold off;

figure(2), imshow(g);
[Image2X,Image2Y] = ginput(8);
hold on;
plot(Image2X,Image2Y);

Euc = ones(8, 1);
Image1 = [Image1X, Image1Y, Euc];
Image2 = [Image2X, Image2Y, Euc];

%Computing the midpoint for all 3 Images

t1 = [mean(Image1X), mean(Image1Y)]; % Midpoint for Image 1
disp("The midpoint for Image 1 is:");
disp(t1);

t2 = [mean(Image2X), mean(Image2Y)]; % Midpoint for Image 2
disp("The midpoint for Image 2 is:");
disp(t2);


% Making the midpoint as origin for the 3 Images and computing Translated
% corrdinate

X1_T = (Image1X - t1(1,1));
Y1_T = (Image1Y - t1(1,2));
Image1Translated = [X1_T, Y1_T]; % Translated points for Image 1


X2_T = (Image2X - t2(1,1));
Y2_T = (Image2Y - t2(1,2));
Image2Translated = [X2_T, Y2_T]; % Translated points for Image 2


% Computing the absolute average of the translated points defined as
% SCALING

s1 = [mean(abs(X1_T)), mean(abs(Y1_T))]; % Scaling for Image 1
disp("The scaling for Image 1 is: ");
disp(s1);


s2 = [mean(abs(X2_T)), mean(abs(Y2_T))]; % Scaling for Image 2
disp("The scaling for Image 2 is: ");
disp(s2);


% Computing the Coordinate Transformation for the 3 Images
% T = S TransMatrix X Midpoint TransMatrix

s1_TransMatrix = [1/(s1(1,1)), 0, 0; 0, 1/(s1(1,2)), 0; 0, 0, 1];
t1_TransMatrix = [1, 0, -1*(t1(1,1)); 0, 1, -1*(t1(1,2)); 0, 0, 1];

Image1Transformation = s1_TransMatrix * t1_TransMatrix; % Coordinate transformation for Image 1
disp("The coordinate transformation for Image 1 is: ");
disp(Image1Transformation);


s2_TransMatrix = [1/(s2(1,1)), 0, 0; 0, 1/(s2(1,2)), 0; 0, 0, 1];
t2_TransMatrix = [1, 0, -1*(t2(1,1)); 0, 1, -1*(t2(1,2)); 0, 0, 1];

Image2Transformation = s2_TransMatrix * t2_TransMatrix; % Coordinate transformation for Image 2
disp("The coordinate transformation for Image 2 is: ");
disp(Image2Transformation);



%Computing the conditioned coordinates of the image

Image1Conditioned = [X1_T * Image1Transformation(1,1), Y1_T * Image1Transformation(2,2)]; % Conditioned Coordinaate points for Image 1
disp("The conditioned coordinate of Image 1 is: ");
disp(Image1Conditioned);

Image2Conditioned = [X2_T * Image2Transformation(1,1), Y2_T * Image2Transformation(2,2)]; % Conditioned Coordinaate points for Image 2
disp("The conditioned coordinate of Image 2 is: ");
disp(Image2Conditioned);


conditioned1 = transpose(Image1Conditioned);
conditioned2 = transpose(Image2Conditioned);

% Construction of the design maxtrix A for Image 1 and 2

A = zeros(1+size(conditioned1, 2), 9);

for i = 1:size(conditioned1,2)
    
    A(i, :) = [conditioned1(1, i)*conditioned2(1, i), conditioned1(2, i)*conditioned2(1, i), conditioned2(1, i), conditioned1(1, i)*conditioned2(2, i), conditioned1(2, i)*conditioned2(2, i), conditioned2(2, i), conditioned1(1, i), conditioned1(2, i), 1];
    
    disp(" The Design Matrix A for Image 1&2 is expressed as:");
    disp(A);
   
end


[U, D, V] = svd(A); %Single value decomposition

disp("The value of D is");
disp(D);

fMatrix = reshape(V(:,end), 3, 3)'; %Computing the Fundamental matrix for projective point for image 1&2

fMatrix = Image2Transformation' * fMatrix * Image1Transformation; %Computing the reverse conditioning

fMatrix = fMatrix(:,:)/fMatrix(end,end);  %Normalizing

fMatrix_P = singularity(fMatrix);


disp("Image 1");
disp(Image1);

disp("Image 2");
disp(Image2);

%Calculating and Drawing the Epipolar Lines for all points in the  Images  

%for Image 1 and 2 the Epipolar Line is Calculated as

for i = 1:size(Image1, 1)
    Image2VectorPoints = Image2(i,:);
    lImage1 = fMatrix_P' * Image2VectorPoints';
    l1(i,:) = lImage1'

    Image1VectorPoints = Image1(i,:);
    lImage2 = fMatrix_P * Image1VectorPoints';
    l2(i,:) = lImage2'
end

disp("Image 1 Epipolar Lines are:");
disp(l1);

disp("Image 2 Epipolar Lines are:");
disp(l2);


%Proceeding to draw the corresponding Lines on the Image

%for Image 1
f = imread('b1.jpg');                                   % Read image f
figure, imshow(f), hold on                     % Show image in new window
for i = 1:size(l1,1)
    ln1 = l1(i,:)';
    hline(ln1);                                    % Draw red line l in image f
end



%for Image 2
f = imread('b2.jpg');                                   % Read image f
figure, imshow(f), hold on                     % Show image in new window
for i = 1:size(l2,1)
    ln2 = l2(i,:)';
    hline(ln2);                                    % Draw red line l in image f
end


%Calculating the geometric error

for i = 1:size(Image2, 1)
    a = Image2(i, 1);
    b = Image2(i, 2);
    

    x = l2(i, 1)';
    y = l2(i, 2)';
    z = l2(i, 3)';

    u = Image1(i, 1);
    v = Image1(i, 2);
    

    numerator = (a*x + b*y + z)^2;
    denominator = a^2 + b^2 + u^2 + v^2;
    geometricError = numerator/denominator;

    disp("The Geometric Error is calculated as:")
    disp(geometricError);
end
end




%Implementing a function to compute signularity constraint

function fMatrix_P = singularity(fMatrix)

if det(fMatrix) == 0

    fMatrix_P = fmatrix;
    disp("The Fundamental matrix is:");
    disp(fMatrix_P);

else

    [U, D, V] = svd(fMatrix);
    
    % replacing the last diagonal element of D with 0

    D(end, end) = 0;
    fMatrix_P = U * D * V';
    disp("The Fundamental matrix is:");
    disp(fMatrix_P);
end 
end






function hline(l, varargin)                       
%        ==================

if abs(l(1)) < abs(l(2))                                  % More horizontal
    xlim = get(get(gcf, 'CurrentAxes'), 'Xlim');  
    x1 = cross(l, [1; 0; -xlim(1)]);
    x2 = cross(l, [1; 0; -xlim(2)]);
else                                                        % More vertical
    ylim = get(get(gcf, 'CurrentAxes'), 'Ylim');
    x1 = cross(l, [0; 1; -ylim(1)]);
    x2 = cross(l, [0; 1; -ylim(2)]);
end
x1 = x1 / x1(3);
x2 = x2 / x2(3);
line([x1(1) x2(1)], [x1(2) x2(2)], varargin{:});
end


  



