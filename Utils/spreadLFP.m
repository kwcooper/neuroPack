function [lfp_,splitMult] = spreadLFP(lfp,splitMult)
% returns lfp with offsets for multi channel plotting


if ~exist('splitMult','var') || isempty(splitMult)
  splitMult = 2.5;
end

% Flip the lfp matrix to keep things uniform
%lfp = flip(lfp);

% compute a porportinate offset
split = splitMult * rms(double(lfp(:)));

[nElecs,tPts,nSets] = size(lfp);
lfp_ = lfp/split;
offsets = repmat([1:nElecs]',1,tPts,nSets);
lfp_ = lfp_ + flipud(offsets); 

%lfp__ = flip(lfp_);
end

