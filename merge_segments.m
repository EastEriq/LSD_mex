function merged=merge_segments(segs,offline,angle,score)
% given a matrix of 7xN segments returned by LSD, try to merge greedily
%  all those which are nearly colinear. The result is dependent on the
%  order of thesegments returned by LSD.
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
    angle  = 1*pi/180;
    score = -Inf;
end

sintol=sin(angle);

% remove segments with low score
segs=segs(:,segs(7,:)>score);

% initialize empty merged
merged=zeros(4,0);

% matched segments are removed from segs till none is left
while ~isempty(segs)

    % pick up the first orphan
    testsegment=segs(1:4,1);

    % find all other segs which are colinear (that will include testsegment
    %  itself): 
    % - check that both extremes of the new segment are within
    %  tolerance on the line of the tested segment,
    % - check that the segments are nearly parallel
    % the check could be made slightly more efficient, computing L3 and s2
    %  only for the segments passing the check on s1
    L1=sqrt((segs(3,:)-segs(1,:)).^2 + (segs(4,:)-segs(2,:)).^2);
    L2=sqrt((segs(1,:)-testsegment(1)).^2 + (segs(2,:)-testsegment(2)).^2);
    % L2(1) is always 0
    L3=sqrt((segs(1,:)-testsegment(3)).^2 + (segs(2,:)-testsegment(4)).^2);
    s1=( (testsegment(3)-testsegment(1))*(testsegment(2)-segs(2,:)) - ...
         (testsegment(4)-testsegment(2))*(testsegment(1)-segs(1,:)) )./(L1(1)*L2);
    s1(1)=0;
    s2=( (testsegment(3)-testsegment(1))*(testsegment(2)-segs(4,:)) - ...
         (testsegment(4)-testsegment(2))*(testsegment(1)-segs(3,:)) )./(L1(1)*L3);
    s3=( (testsegment(3)-testsegment(1))*(segs(4,:)-segs(2,:)) - ...
         (testsegment(4)-testsegment(2))*(segs(3,:)-segs(1,:)) )./(L1(1)*L1);
    candidates=find(abs(s1)<sintol & abs(s2)<sintol & abs(s3)<sintol);

    % project the extremes of each candidate segment on testsegment
    p1 = segs(1,candidates)*(testsegment(3)-testsegment(1)) + ...
         segs(2,candidates)*(testsegment(4)-testsegment(2));
    p2 = segs(3,candidates)*(testsegment(3)-testsegment(1)) + ...
         segs(4,candidates)*(testsegment(4)-testsegment(2));

    % reorient candidate segments so that p1<p2 always
    toreverse=(p1>p2);
    temp=segs(1,candidates(toreverse));
    segs(1,candidates(toreverse))=segs(3,candidates(toreverse));
    segs(3,candidates(toreverse))=temp;
    temp=segs(2,candidates(toreverse));
    segs(2,candidates(toreverse))=segs(4,candidates(toreverse));
    segs(4,candidates(toreverse))=temp;
    
    % take the two extremal points as extremes of the merged segment
    % This implies overly relying on their plausibility. Maybe some kind of
    %  fit which takes all the intermediate points into account would be
    %  better?
    [~,i1]=min(p1);
    [~,i2]=max(p2);
    X1=segs(1:2,candidates(i1));
    X2=segs(3:4,candidates(i2));

    merged=[merged,[X1;X2]];

    % remove from segs all the candidates used
    segs(:,candidates)=[];
end