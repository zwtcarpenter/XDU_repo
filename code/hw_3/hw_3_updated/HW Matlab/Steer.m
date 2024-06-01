
function X_new=Steer(X_rand,X_near,StepSize)
theta = atan2(X_rand(1)-X_near(1),X_rand(2)-X_near(2));  % direction to extend sample to produce new node
X_new = X_near(1:2) + StepSize * [sin(theta)  cos(theta)];