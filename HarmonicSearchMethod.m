clc;
clear;
close all;
%% Problem Definition
CostFunction=@cost;        % Cost Function
nVar=1;                    % Number of Decision Variables
xVarSize=[1 nVar];         % Decision Variables Matrix Size
xVarMin= 5492;              % Decision Variables Lower Bound
xVarMax= 6307;               % Decision Variables Upper Bound
%% Harmony Search Parameters
MaxIt=500;                 % Maximum Number of Iterations
HMS=50;                    % Harmony Memory Size
nNew=100;                  % Number of New Harmonies
HMCR=0.9;                  % Harmony Memory Consideration Rate
PAR=0.1;                   % Pitch Adjustment Rate
FW=0.02*(xVarMax-xVarMin); % Fret Width (Bandwidth)
%% Initialization

% Empty Harmony Structure
empty_harmony.Position=[];
empty_harmony.Cost=[];

% Initialize Harmony Memory
HM=repmat(empty_harmony,HMS,1);

% Create Initial Harmonies
for i=1:HMS
    rand
    HM(i).Position=unifrnd(xVarMin,xVarMax,xVarSize);
    HM(i).Cost=CostFunction(HM(i).Position);
end

% Sort Harmony Memory
[~, SortOrder]=sort([HM.Cost]);

% Update Best Solution Ever Found
BestSol=SortOrder(1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

%% Harmony Search Main Loop
for it=1:MaxIt
    
    % Initialize Array for New Harmonies
    NEW=repmat(empty_harmony,nNew,1);
    
    % Create New Harmonies
    for k=1:nNew
        
        % Create New Harmony Position
        NEW(k).Position=unifrnd(xVarMin,xVarMax,xVarSize);

        for j=1:nVar
            if rand<=HMCR
                % Use Harmony Memory
                i=randi([1 HMS]);
                NEW(k).Position(j)=HM(i).Position(j);
            end
            
            % Pitch Adjustment
            if rand<=PAR
                DELTA=FW*randn();    
                NEW(k).Position(j)=NEW(k).Position(j)+DELTA;
            end
        
        end
        
        % Apply Variable Limits
        NEW(k).Position=max(NEW(k).Position,xVarMin);
        NEW(k).Position=min(NEW(k).Position,xVarMax);
        
        % Evaluation
        NEW(k).Cost=CostFunction(NEW(k).Position);
        
    end
    
    % Merge Harmony Memory and New Harmonies
    HM=[HM 
        NEW]; 
    
    % Sort Harmony Memory
    [~, SortOrder]=sort([HM.Cost]);
    HM=HM(SortOrder);
    
    % Truncate Extra Harmonies
    HM=HM(1:HMS);
    
    % Update Best Solution Ever Found
    BestSol=HM(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        
end
%% Results
figure;
plot(BestCost,'LineWidth',2);
%semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Optimal Cost');
grid on;
