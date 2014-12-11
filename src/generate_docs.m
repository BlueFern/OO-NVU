%% Script for generating documentation
opts.outputDir = '../doc/html';
opts.showCode = false;
publish('astrocyte_doc.m', opts);
publish('smcec_doc.m', opts);
publish('wallmechanics_doc.m', opts);
opts.showCode = true;
publish('nvu_script.m', opts);