# CoreForge — Örnek Çalışmalar (Adım Adım Kılavuzlar)

Şartnamedeki "örnekler ve kılavuzlar" şartını karşılayan altı uygulamalı
çalışma. Her biri 2–10 dakika sürer; beklenen sayısal sonuçlar verilir,
böylece kurulumunuzu da sınamış olursunuz.

---
## Örnek 1 · SMR korunun 3-D analizi (ana kullanım senaryosu)

1. Uygulama SMR koruyla açılır (📚'dan **SMR-class PWR** de yükleyebilirsiniz).
2. Kenar çubuğunda **Solver model → 3-D** seçin: kor tek zon olarak
   extrude edilir (10 katman × 20 cm), alt/üst BC eklenir.
3. 🧱 Core builder → zon yöneticisi: aktif zonu seçin, **➕ below** ile
   alta kopya ekleyin, yeni zonda **▩ fill = 5** (ağır reflektör),
   layers = 1. Aynısını **➕ above** ile üste yapın.
4. ⚡ Solve. Beklenen: k_eff ≈ 1.04 mertebesi; **eksenel güç profili**
   kosinüs benzeri, F_z ≈ 1.4–1.6; katman gezgini ile z-dilimlerini
   inceleyin.
5. **🔌 Operating point**: P=160 MW girin → q‴ ve termal akı
   (~10¹³–10¹⁴ n/cm²·s bandında) mutlak birimlerde raporlanır.
6. 💾 Save project — bu sizin 3-D SMR tasarım dosyanızdır.

## Örnek 2 · IAEA-3D benchmark'ı ve S-eğrisi

1. 📚 → **IAEA-3D PWR** → Load → ⚡ Solve (varsayılan mesh, <1 s).
   Beklenen: k=1.028741, referans farkı ≈ −27 pcm (kenar çubuğundan
   div=4, divz=4 yaparsanız ≈ −13 pcm; süre birkaç saniye).
2. 🔧 Physics tools → **2 · Rod insertion sweep** → 6 derinlik → çalıştır.
   Beklenen: toplam çubuk-5 değeri ≈ −980 pcm, S-şekilli eğri.
3. Yorum: eğrinin orta dikliği eksenel akı tepesiyle çakışır — yalnız
   3-D'nin gösterebildiği fizik.

## Örnek 3 · Yakıt tasarımı ve kritik boron

1. 📚 → **Designer PWR** → Load → ⚡ Solve (k≈1.105 @1000 ppm).
2. 🔧 → **3 · Critical boron search** → hedef k=1.0 → çalıştır.
   Beklenen: **≈ 2 236 ppm**, ~16 koşu.
3. 🧬'de %4.5'lik yakıt üretin, **Add to materials**; 🧱'de merkeze
   boyayın; tekrar Solve → k artışını ve güç haritası kaymasını izleyin.

## Örnek 4 · Yakıt çevrimi: letdown ve otomatik çevrim sonu

1. Designer PWR yüklüyken 🔥 Burnup → **auto end-of-cycle** işaretli,
   boron = *critical letdown* → **Run fuel cycle** (~30–60 s).
2. Beklenen: letdown 1 767 → 0 ppm (her adım k=1.0000), çevrim
   **30.0 MWd/kgU = 789 EFPD**, "reactivity-limited" rozeti; BOC xenon
   değeri ≈ −3 570 pcm; EOC yanma haritası 5–26 MWd/kgU bandında.
3. cycle_history.csv'yi indirin — rapor eklerinize hazırdır.

## Örnek 5 · Benchmark yakıtını yakılabilir yapmak (ters çözücü)

1. 📚 → **IAEA-2D** → Load. 🧬 → **🔁 Equivalent fuel** → "Fuel 1" →
   **Fit**. Beklenen: e=2.886 w/o, 259 ppm; νΣf₂/Σa₂ rezidüsü ~0.
2. **Attach** → malzeme "(eq-fuel)" olur. Artık 🔥 Burnup bu 1976
   benchmark koru üzerinde çalışır (yaklaşık eşdeğer, etiketli).

## Örnek 6 · Transient: çubuk fırlatması ve scram

1. ⏱ Transient → P₀=160 MW, **step**, ρ=+300 pcm, **Doppler feedback
   açık** (α_D=−2.5 pcm/K) → Run.
   Beklenen: güç ~**316 MW tepe** yapıp öz-sınırlanır; yakıt sıcaklığı
   580→~685 K; net ρ → +37 pcm.
2. Aynı ayarlarla feedback'i KAPATIN, ρ=+100 pcm, 120 s: asimptotik
   periyot ≈ **54.9 s** — panel inhour değeriyle karşılaştırır.
3. **scram**, ρ=−4 000 pcm, gecikme 1 s: 5 s'te ~%7 güce iniş
   (bozunum ısısı modellenmez — kinetik seviye).

---
## Örnek 7 · Isıl-hidrolik, iyot çukuru ve çok-çevrim

1. SMR yüklüyken ⚡ Solve, sonra 🌡 **Thermal-hydraulics** → varsayılan
   girdilerle **Run channel analysis**.
   Beklenen: çıkış sıcaklığı ~269/273 °C (ort/sıcak), yakıt merkez
   ~820 °C, doyma marjı ~+61 K, ΔP ~15 kPa; eksenel sıcaklık
   profillerini inceleyin. Enerji kapanışı metriği kor gücüne eşittir.
2. ⏱ Transient → **☁️ Xenon transient** → %0 güç (kapanma) → Run.
   Beklenen: iyot çukuru ~7-8 saatte, denge değerinin ~×1.5'i;
   toparlanma ~22 saat. "Restart penceresi" kavramını grafik üzerinden
   okuyun.
3. 📚 → Designer PWR → 🔥 Burnup → auto + letdown → Run (~1 dk).
   Bitince **🔁 Multi-cycle** genişleticisi → fresh fraction 0.34 →
   **Run next cycle**. Beklenen: 2. çevrim BOC boronu 1. çevrimden
   düşük, çevrim uzunluğu kısalır; tablo çevrim çevrim birikir —
   denge çevrimine yaklaşımı izleyin.

---
Her örneğin sonunda ⚡ sekmesindeki **📄 HTML report** ile tek dosyalık
arşiv alınabilir; 💾 proje dosyası çalışmayı aynen geri yükler.
