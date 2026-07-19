# CoreForge — Bağımsız Doğrulama Protokolü

Bu belge, CoreForge'un sonuçlarını **koddan bağımsız** bir denetçiye (bir
insan gözden geçiren veya ayrı bir dil modeli oturumu) yeniden ürettirmek,
sonra kodu **düşmanca** sınamak için hazır komut istemleri içerir. Amaç,
"kendi kendini onaylayan" doğrulamadan kaçınmaktır: sayıları üreten taraf ile
onları sorgulayan taraf farklı olmalıdır.

**Temel pratik:** Bu turlarda bulunan HER somut hata/eksik, `verify.py`'ye
**kalıcı bir kontrole** dönüştürülür. Doğrulama tek seferlik bir olay değil,
regresyon kapısıdır.

> Denetçiye kod tabanının tamamına salt-okunur erişim verin. Denetçinin
> `python3 verify.py` (derleyici varsa) veya `python3 verify.py --no-engine`
> (derleyici yoksa, saf-Python altküme) çalıştırabildiğinden emin olun.

---

## Tur 1 — Çalıştırmalı yeniden üretim (reproduction)

Amaç: dokümandaki her sayının gerçekten bu koddan çıktığını bağımsızca
doğrulamak.

**Hazır komut istemi:**
```
Sana bir nötron difüzyon reaktör analiz kodu (CoreForge) verildi. Hiçbir
sonuca güvenme; hepsini kendin çalıştırarak doğrula.

1. Motoru derle:
     ifx -O3 -qopenmp solver/coreforge.f90 -o solver/coreforge
     (ifx yoksa: gfortran -O3 -fopenmp ...  ; hiç derleyici yoksa 3. adıma geç)
2. `python3 verify.py --fine` çalıştır. Tüm kontroller PASS mı?
3. Derleyici yoksa `python3 verify.py --no-engine` çalıştır.
4. docs/DOGRULAMA_RAPORU.md ve README.md'deki k_eff ve pcm-fark tablolarını
   al; verify.py çıktısındaki gerçek sayılarla SATIR SATIR karşılaştır.
   Uyuşmayan tek bir hücre var mı? Varsa dosya:satır ile bildir.
5. IAEA-2D için mesh yakınsamasını kendin kontrol et: div=5, 8, 10 ile
   çalıştır; k_eff monoton mu ve referansa (1.02959) yakınsıyor mu?

Çıktı: her sayı için "yeniden üretildi / sapma X pcm / üretilemedi" tablosu.
```

**Beklenen:** Tüm kontroller PASS; doküman sayıları koşu sayılarıyla ≤ birkaç
pcm uyumlu (mesh/optimizasyon kaynaklı küçük farklar kabul edilebilir, işaret
ve büyüklük mantıklı olmalı).

---

## Tur 2 — Düşmanca sınama (adversarial)

Amaç: kodu kırmaya, gizli varsayımları ve sessiz hataları ortaya çıkarmaya
çalışmak. Denetçi "bunu nasıl yanlış sonuca zorlarım?" diye düşünmeli.

**Hazır komut istemi:**
```
Bu kodun YANLIŞ olduğunu kanıtlamaya çalış. Şu saldırı vektörlerini dene ve
her biri için somut bir girdi + gözlemlenen davranış bildir:

- Fizik değişmezlikleri: mesh inceltme (div/divz artışı) k_eff'i referans
  toleransının ötesinde KAYDIRIYOR mu? (Kaydırıyorsa ayrıklaştırma hatası
  vardır.) QA parmakizi mesh değişiminde SABİT, ama pitch/kesit/BC
  değişiminde DEĞİŞİYOR mu? (presets.physics_fingerprint)
- Sınır koşulları: tümü yansıtıcı (reflective) BC'de sonsuz-ortam k_inf'i
  analitik nuSf/Sa değerini veriyor mu? Vakum BC gevşetilince k düşüyor mu?
- Simetri: simetrik bir çekirdekte akı/güç haritası simetrik mi? Köşe/kenar
  hücre indekslemesinde kayma var mı? (geçmiş hata: burnup pmap köşe hücresi)
- Korunum: T-H enerji dengesi kapanıyor mu (Q = ṁ·cp·ΔT)? Akı kalibrasyonu
  belirtilen gücü geri veriyor mu (S = q_avg/κ, ν̄ İÇERMEMELİ)?
- Birim tuzakları: ksenon σφ birimi (0.053 vs 5.3e-5), fz>π/2 kosinüs clamp
  bayrağı sessizce mi tetikleniyor yoksa raporlanıyor mu?
- Sayısal uçlar: 1 grup vs çok grup, NZ=1 (2B) ile tek-katman extrude 3B
  bit-düzeyinde aynı mı? Yakınsamayan bir vaka "converged" diyor mu?

Her bulguyu: {dosya:satır, tetikleyen girdi, beklenen, gözlemlenen, önem}
biçiminde ver.
```

**Beklenen:** Kritik bir hata bulunmaması idealdir; ama bulunursa bu bir
başarıdır — düzeltilir ve Tur 3'e geçilir.

---

## Tur 3 — Referans kontrol (independent-source cross-check)

Amaç: kodun ürettiği değerleri, kodla **hiçbir ortak kökeni olmayan** dış
kaynaklarla karşılaştırmak.

**Hazır komut istemi:**
```
CoreForge'un sonuçlarını bağımsız kaynaklarla çapraz-doğrula:

1. IAEA-2D/3D benchmark: yayınlanmış referans (Argonne ANL-7416 / IAEA
   benchmark seti) k_eff = 1.02959 (2B). Kodun sonucu bu değere kaç pcm
   uzakta ve fark kaynağı (ayrıklaştırma mı, model mi) açıklanabiliyor mu?
2. OECD/NEA C5G7 benchmark: referans k_eff = 1.18655. Koddaki homojenize
   difüzyon yaklaşımının bu değerden sapması (~+31 pcm Richardson ekstrapo-
   lasyonlu) fiziksel olarak makul bir MODEL hatası mı, yoksa uygulama hatası
   mı? c5g7_data.py kesitleri referansla tutarlı mı?
3. Yakıt tasarımcısının (xslib) ürettiği 2-grup kesitleri, kodla ilgisiz
   IAEA-2D yakıt sabitleriyle (D1≈1.5, D2≈0.4, Sa1≈0.010, Sa2≈0.080,
   S12≈0.020) kıyasla — mühendislik mertebesinde tutuyor mu?
4. Nokta kinetiği: +100/−100 pcm adım için periyodu analitik inhour
   denklemiyle bağımsız hesapla; kodun periyoduyla %1 içinde mi?

Çıktı: her kaynak için {kod değeri, dış referans, fark, "model hatası /
uygulama hatası / kabul" yargısı}.
```

**Beklenen:** Sapmalar ya analitik olarak açıklanabilir (C5G7 model sınırı
gibi) ya da tolerans içinde olmalı.

---

## Bulguları kalıcılaştırma (regresyon kapısı)

Herhangi bir turda **somut bir hata** bulunursa:

1. Kök nedeni düzelt.
2. `verify.py`'ye, o hatayı bir daha yakalayacak **yeni bir kontrol** ekle
   (mevcut kontrollerin biçimini izle: bir `*_check()` fonksiyonu, `PASS/FAIL`
   basar, `main()` içinde hem tam hem `--no-engine` uygun yola bağlanır).
3. Kontrol sayısını güncelle: `README.md` ve `docs/DOGRULAMA_RAPORU.md`.
4. `python3 verify.py --fine` yeniden yeşil olana kadar tekrarla.

Bu döngü projenin temel kalite pratiğidir: her bulunan hata, kalıcı bir
otomatik kontrole dönüşür.

---

## Örnek: geçmişte bu yolla eklenmiş kontroller

- **Sürüm-lockstep** (`version_consistency_check`): motor/report/app sürüm
  dizeleri birbirinden kaydığında sert FAIL. (Yayın öncesi motor 8.2 iken
  arayüz 8.3 kayması bu kontrolle yakalandı.)
- **QA parmakizi** (`integrity_guard_check`): mesh inceltme referansı korur,
  fizik değişimi (pitch/kesit/BC) referansı geçersiz kılar.
- **Akı kalibrasyonu** (`flux_calibration_check`): enerji kapanışı, raporda
  ν̄ fazlalığı hatasına karşı.
- **Denge yakıt güvenilirlik bayrağı** (`eqfuel_reliability_check`): güçlü
  soğurucu (kontrol çubuğu) malzemesinde uyarı, düz yakıtta sessiz.
