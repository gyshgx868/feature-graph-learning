%% Create currupt version of the model
% Takes a shape filename (assumed to be of the name XXX_GT.mat), and
% creates noisy versions of it.
% options is a structure with the following optional fields:
% 'GaussianNoise' - a list of standard deviations for Gaussian additive corrupted versions of the surface. The standard deviation is expressed as a fraction of the surface' diameter
% 'ShotNoiseProbability','ShotNoiseNumber','ShotNoiseMagnitude' - parameter lists for shot noise (i.e point are replaced by arbitrary 3D values. Parameter lists should be of equal size
% 'DepthShotNoiseProbability','DepthShotNoiseMagnitude','DepthShotNoiseNumber' - for depth shot noise (i.e point are perturbed in the direction of the ray leading from the range scanner to the object
% 'CameraCenter' - assumed camera location, for depth noise

function X = addnoise(X_gt, shapename, options)
if (exist('options','var')==0)
    options=[];
end
% defaults=struct('GaussianNoise',[0.00125 0.0025 0.005 0.01 ],...
%     'ShotNoiseProbability',[0 0],...
%     'ShotNoiseNumber',[10 100],...
%     'ShotNoiseMagnitude',[1 1]*0.2,...
%     'DepthShotNoiseProbability',0,...
%     'DepthShotNoiseMagnitude',0.2,...
%     'DepthShotNoiseNumber',100,...
%     'CameraCenter',[ 0 0 -6]); % optical center for GIP camera
% options=incorporate_defaults(options,defaults);

SAMPLING_SET=200;
srf=struct('X',X_gt(:,1),'Y',X_gt(:,2),'Z',X_gt(:,3));
% estimate diameter
ifps=fps_euc(srf,SAMPLING_SET);
Dfps=pdist2(X_gt(ifps,:),X_gt(ifps,:));
diam=sqrt(max(Dfps(:)))

% % Create surface with holes
% for i=1:length(options.DepthShotNoiseProbability)
%     prb=options.DepthShotNoiseProbability(i);
%     X=X_gt;
%     if (prb>0)
%         idx=find(rand(size(X_gt(:,1)))<prb);
%         prb_name=num2str(prb);
%     else
%         idx=randperm(size(X_gt,1),options.DepthShotNoiseNumber(i));
%         prb_name=num2str(options.DepthShotNoiseNumber(i));
%     end
%     dr=(X_gt(idx,:)-repmat(options.CameraCenter(:)',[numel(idx),1]));
%     ndr=sqrt(sum(dr.^2,2)+1e-6);
%     X(idx,:)=X_gt(idx,:)+(repmat(randn(size(X_gt(idx,1))),[1 3]).*bsxfun(@rdivide,dr,ndr)*diam*options.DepthShotNoiseMagnitude(i));
%     mat_filename=[shapename,'_depth_shot_noise_',prb_name,'.mat'];
%     save(mat_filename,'X');
%     ply_filename=[shapename,'_depth_shot_noise_',prb_name,'.ply'];
%     write_ply_only_points(X,ply_filename);
% end

% Create a Gaussian noised version of the surface
%for i=1:length(options.GaussianNoise)
    sig=diam*options.GaussianNoise
    X=X_gt+randn(size(X_gt))*sig;    
    %mat_filename=[shapename,'_gaussian_noise_',num2str(sig),'.mat'];
    %save(mat_filename,'X');    
    ply_filename=[shapename,'_gaussian_noise_',num2str(options.GaussianNoise)];        
    %ply_filename=[shapename,'_gaussian_noise_',num2str(sig),'.ply'];
    %#####################################
    pcwrite(pointCloud(X), ['models/noise/' ply_filename '.ply']);
    %####################################
%end

% % Create a 3D shot noise version of the surface
% for i=1:length(options.ShotNoiseProbability)
%     prb=options.ShotNoiseProbability(i);
%     X=X_gt;
%     if (prb>0)
%         idx=find(rand(size(X_gt(:,1)))<prb);
%         prb_name=num2str(prb);
%     else
%         idx=randperm(size(X_gt,1),options.ShotNoiseNumber(i));
%         prb_name=num2str(options.ShotNoiseNumber(i));
%     end
%     X(idx,:)=X_gt(idx,:)+randn(size(X_gt(idx,:)))*diam*options.ShotNoiseMagnitude(i);
%     mat_filename=[shapename,'_shot_noise_',prb_name,'.mat'];
%     save(mat_filename,'X');
%     ply_filename=[shapename,'_shot_noise_',prb_name,'.ply'];
%     write_ply_only_points(X,ply_filename);
% end

end
