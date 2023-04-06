function [s, w]=runFslCommand(rap,fslcmd,ENV)

% Setup
fslsetup = deblank(rap.directoryconventions.fslsetup);
if not(isempty(fslsetup))
    if ~startsWith(fslsetup,'. '), fslsetup = ['. ' fslsetup]; end
    if ~endsWith(fslsetup,';'), fslsetup=[fslsetup ';']; end
end

if nargin < 3, ENV = {}; end

global reproacache
if ~reproacache.isKey('toolbox.spm'), logging.error('SPM is not found'); end
SPMtool = reproacache('toolbox.spm');

ENV = vertcat(ENV,{...
    'FSLOUTPUTTYPE', rap.directoryconventions.fsloutputtype; ...
    'MATLABPATH', SPMtool.toolPath;...
    });

[s, w] = shell(fslcmd,'shellprefix',fslsetup,'environment',ENV);

% Display error if there was one
if s, logging.error('runFslCommand:Error running %s\n%s',cmd,w); end