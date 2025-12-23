function merged=merge_segments(segs,offline,angle,score)
% given a matrix of 7xN segments returned by LSD, try to merge greedily
%  all those which are nearly colinear
% LSD at high image resolution (no subsampling) is known to return many
%  short segments instead of long lines; moreover, in our typical
%  astronomical images we can safely assume that streaks are long and not
%  intermittent
% LSD is essentially local in identifying segment candidates, and weak
%  streaks are also intermittently detected because of noise
%
% Inputs: segs     - only segs([1:4,7},:) are used
%         offline  - offset threshold of a segment point from the searched
%                    line. If the start point of the second segment is less
%                    than that, the segment is merged in
%         angle    - angular tolerance to group together segments
%         score    - threshold on -log(NFA) (segs(7,:)). Input segments
%                    below that score are neglected
% Output: merged   - an array of 4xN segment coordinates (no attempt is
%                    made to return a merged value of the corresponding
%                    values of segs(5:7,:) )
arguments
    segs double
    offline = 2;
    angle  = 2*180/pi;
    score = -Inf;
end

% remove segments with low score
segs=segs(:,segs(7,:)>score);

% initialize empty merged
merged=zeros(4,0);

% matched segments are removed from segs till none is left
while ~isempty(segs)
    % pick up the first orphan
    testsegment=segs(1:4,1);
    % find all other segs which are colinear (that will include testsegment
    %  itself)
    L=sqrt((segs(3,:)-segs(1,:)).^2 + (segs(4,:)-segs(2,:)).^2);
    s1=(testsegment(3)-testsegment(1))*(segs(4,:)-segs(2,:)) - ...
       (testsegment(3)-testsegment(1))*(segs(4,:)-segs(2,:))
    candidates=find();
    % project the extremes of each candidate segment on testsegment
    
    % sort all the candidates according to the projected distance on the
    %  line
    
    % take the two extremal points as extremes of the merged segment
    merged=[merged,];
    % remove from segs all the candidates considered
    segs(:,candidates)=[];
end