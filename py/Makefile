BARDIR=../bar

bar: $(BARDIR)/rphi2d.ps $(BARDIR)/rphi1d.ps $(BARDIR)/detecta.ps \
	$(BARDIR)/detectb.ps $(BARDIR)/detectc.ps $(BARDIR)/detectd.ps \
	$(BARDIR)/rphi1d2.ps

$(BARDIR)/rphi2d.ps: bar_figures.py calc_veldist_2d.py
	$(PYTHON) bar_figures.py $@ --saveDir='../bar/2d'

$(BARDIR)/rphi1d.ps: bar_figures.py calc_veldist_1d.py
	$(PYTHON) bar_figures.py $@ --vlos --saveDir="../bar/1d"

$(BARDIR)/rphi1d2.ps: bar_figures.py calc_veldist_1d.py
	$(PYTHON) bar_figures.py $@ --vlos --saveDir="../bar/1d2" --otherside

$(BARDIR)/detecta.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@

$(BARDIR)/detectb.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@ --skipCenter=0.15

$(BARDIR)/detectc.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@ --convolve=0.1

$(BARDIR)/detectd.ps: bar_detectability.py
	$(PYTHON) bar_detectability.py $@ --skipCenter=0.15  --convolve=0.1

